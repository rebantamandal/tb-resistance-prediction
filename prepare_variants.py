"""Convert a long-format M. tuberculosis variant call table into the isolate-level
feature matrix that `amr.py` consumes.

Input  : all_variants.csv with columns SAMPLE,CHROM,POS,REF,ALT (one row per
         variant call per isolate; ~15M rows, ~9.8k isolates).
Output : a wide binary matrix, one row per isolate, one 0/1 column per retained
         variant, plus a feature dictionary annotating each column with the
         H37Rv gene it falls in and the drugs that gene is associated with.

This step performs no label inference. Resistance phenotypes must be supplied
separately (--phenotypes); without them the matrix is written unlabelled and is
usable only for prediction, not training. Deriving a resistance label from the
same variants that serve as predictors would be circular and is not supported.

Run `python prepare_variants.py --help` for options.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

CHUNK_ROWS = 2_000_000
VARIANT_COLUMNS = ["SAMPLE", "CHROM", "POS", "REF", "ALT"]
DEFAULT_ORGANISM = "Mycobacterium tuberculosis"
MAX_ALLELE_IN_NAME = 12


# --------------------------------------------------------------------------
# Locus annotation
# --------------------------------------------------------------------------

class LocusIndex:
    """Maps an H37Rv coordinate to a gene, including an upstream/downstream margin.

    The margin exists because several important resistance determinants are
    promoter variants outside the coding sequence (the canonical example being
    the fabG1/inhA promoter for isoniazid and ethionamide). Overlapping windows
    are resolved by taking the first match in file order.
    """

    def __init__(self, records: list[dict[str, object]], margin: int) -> None:
        self.margin = margin
        self.records = records
        self._starts = np.array([int(r["start"]) - margin for r in records], dtype=np.int64)
        self._ends = np.array([int(r["end"]) + margin for r in records], dtype=np.int64)
        self._genes = np.array([str(r["gene"]) for r in records], dtype=object)

    @classmethod
    def from_csv(cls, path: Path, margin: int) -> "LocusIndex":
        frame = pd.read_csv(path, dtype={"gene": str, "drugs": str, "note": str})
        required = {"gene", "start", "end", "drugs"}
        missing = required.difference(frame.columns)
        if missing:
            raise SystemExit(f"Locus table {path} is missing columns: {sorted(missing)}")
        if frame.empty:
            raise SystemExit(f"Locus table {path} contains no loci.")
        records = frame.to_dict("records")
        for record in records:
            if int(record["start"]) > int(record["end"]):
                raise SystemExit(f"Locus {record['gene']} has start greater than end.")
        return cls(records, margin)

    def gene_of(self, positions: np.ndarray) -> np.ndarray:
        """Vectorised lookup. Returns an object array of gene names, '' when unannotated."""
        out = np.empty(positions.shape[0], dtype=object)
        out[:] = ""
        unassigned = np.ones(positions.shape[0], dtype=bool)
        for index in range(self._starts.shape[0]):
            if not unassigned.any():
                break
            hit = unassigned & (positions >= self._starts[index]) & (positions <= self._ends[index])
            out[hit] = self._genes[index]
            unassigned &= ~hit
        return out

    def drugs_by_gene(self) -> dict[str, str]:
        return {str(r["gene"]): str(r.get("drugs", "")) for r in self.records}

    def notes_by_gene(self) -> dict[str, str]:
        return {str(r["gene"]): str(r.get("note", "")) for r in self.records}


# --------------------------------------------------------------------------
# Variant keys and column names
# --------------------------------------------------------------------------

def variant_keys(frame: pd.DataFrame) -> pd.Series:
    """Stable textual key for a variant: '<pos>_<ref>_<alt>'.

    CHROM is deliberately excluded from the key; the caller verifies that the
    file contains exactly one reference contig, so including it would only
    lengthen every column name.
    """
    return (frame["POS"].astype(str) + "_" + frame["REF"].astype(str)
            + "_" + frame["ALT"].astype(str))


def _shorten(allele: str) -> str:
    """Keep long indel alleles out of column names without losing uniqueness.

    The digest must come from hashlib, not the built-in hash(): Python salts
    string hashing per process, so the built-in would give the same allele a
    different column name on every run. That silently produced mutually
    incompatible feature panels across drugs, which is exactly the kind of
    irreproducibility that makes saved models unusable together.
    """
    allele = re.sub(r"[^A-Za-z0-9]", "", str(allele)).upper()
    if len(allele) <= MAX_ALLELE_IN_NAME:
        return allele or "NONE"
    digest = hashlib.blake2b(allele.encode("ascii"), digest_size=3).hexdigest()
    return f"{allele[:MAX_ALLELE_IN_NAME]}x{digest}"


def column_name(gene: str, key: str) -> str:
    position, reference, alternate = key.split("_", 2)
    prefix = re.sub(r"[^A-Za-z0-9]", "", gene) if gene else "intergenic"
    return f"g_{prefix}_{position}_{_shorten(reference)}_{_shorten(alternate)}"


# --------------------------------------------------------------------------
# Pass 1: prevalence
# --------------------------------------------------------------------------

def scan_prevalence(variants_path: Path, loci: LocusIndex | None,
                    restrict_to_loci: bool, progress: bool) -> tuple[Counter, dict[str, str], list[str]]:
    """Count how many distinct isolates carry each variant.

    Returns the prevalence counter, the gene assignment per variant key, and the
    sorted list of isolate identifiers seen in the file.
    """
    prevalence: Counter = Counter()
    gene_of_key: dict[str, str] = {}
    samples: set[str] = set()
    contigs: set[str] = set()
    rows_read = 0

    reader = pd.read_csv(variants_path, chunksize=CHUNK_ROWS, dtype=str,
                         usecols=VARIANT_COLUMNS, keep_default_na=False, na_filter=False)
    for chunk in reader:
        rows_read += len(chunk)
        contigs.update(chunk["CHROM"].unique().tolist())
        samples.update(chunk["SAMPLE"].unique().tolist())

        if restrict_to_loci and loci is not None:
            positions = pd.to_numeric(chunk["POS"], errors="coerce").to_numpy(dtype="float64")
            if np.isnan(positions).any():
                raise SystemExit("POS contains non-numeric values; clean the variant file first.")
            genes = loci.gene_of(positions.astype(np.int64))
            keep = genes != ""
            chunk = chunk.loc[keep]
            genes = genes[keep]
        elif loci is not None:
            positions = pd.to_numeric(chunk["POS"], errors="coerce").to_numpy(dtype="float64")
            genes = loci.gene_of(positions.astype(np.int64))
        else:
            genes = np.array([""] * len(chunk), dtype=object)

        if chunk.empty:
            continue

        keys = variant_keys(chunk)
        # A single isolate can legitimately appear twice for the same variant if
        # the source merged call sets; count isolates, not call records.
        pairs = pd.DataFrame({"SAMPLE": chunk["SAMPLE"].to_numpy(), "key": keys.to_numpy()})
        pairs = pairs.drop_duplicates()
        prevalence.update(pairs["key"].tolist())

        # Assign genes once per distinct key rather than once per call row.
        distinct = pd.DataFrame({"key": keys.to_numpy(), "gene": genes}).drop_duplicates("key")
        for key, gene in zip(distinct["key"].to_numpy(), distinct["gene"].to_numpy()):
            if key not in gene_of_key:
                gene_of_key[key] = gene

        if progress:
            print(f"  pass 1: {rows_read:,} rows, {len(prevalence):,} distinct variants",
                  file=sys.stderr, flush=True)

    if len(contigs) != 1:
        raise SystemExit(f"Expected exactly one reference contig, found {sorted(contigs)}. "
                         "Variant keys would be ambiguous across contigs.")
    if not samples:
        raise SystemExit("No isolates found in the variant file.")
    return prevalence, gene_of_key, sorted(samples)


# --------------------------------------------------------------------------
# Pass 2: matrix
# --------------------------------------------------------------------------

def build_matrix(variants_path: Path, sample_index: dict[str, int],
                 feature_index: dict[str, int], progress: bool) -> np.ndarray:
    matrix = np.zeros((len(sample_index), len(feature_index)), dtype=np.uint8)
    rows_read = 0

    reader = pd.read_csv(variants_path, chunksize=CHUNK_ROWS, dtype=str,
                         usecols=VARIANT_COLUMNS, keep_default_na=False, na_filter=False)
    for chunk in reader:
        rows_read += len(chunk)
        keys = variant_keys(chunk)
        columns = keys.map(feature_index)
        keep = columns.notna()
        if keep.any():
            rows = chunk.loc[keep, "SAMPLE"].map(sample_index).to_numpy(dtype=np.int64)
            cols = columns[keep].to_numpy(dtype=np.int64)
            matrix[rows, cols] = 1
        if progress:
            print(f"  pass 2: {rows_read:,} rows", file=sys.stderr, flush=True)
    return matrix


# --------------------------------------------------------------------------
# Phenotypes
# --------------------------------------------------------------------------

POSITIVE_LABELS = {"1", "r", "resistant"}
NEGATIVE_LABELS = {"0", "s", "susceptible", "sensitive"}


def load_phenotypes(path: Path, sample_column: str, drug_column: str) -> pd.DataFrame:
    """Read a per-isolate phenotype table and normalise one drug column to 0/1.

    Intermediate, unknown and missing outcomes are dropped rather than recoded,
    matching the target policy `amr.py` enforces: an ambiguous DST result is not
    silently treated as susceptible.
    """
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
    for needed in (sample_column, drug_column):
        if needed not in frame.columns:
            raise SystemExit(f"Phenotype file has no column '{needed}'. "
                             f"Available: {list(frame.columns)}")
    frame = frame[[sample_column, drug_column]].copy()
    frame.columns = ["isolate_id", "raw_label"]

    normalised = frame["raw_label"].astype(str).str.strip().str.lower()
    resistant = pd.Series(np.nan, index=frame.index, dtype="float64")
    resistant[normalised.isin(POSITIVE_LABELS)] = 1.0
    resistant[normalised.isin(NEGATIVE_LABELS)] = 0.0

    dropped = int(resistant.isna().sum())
    frame["resistant"] = resistant
    frame = frame.dropna(subset=["resistant"])
    frame["resistant"] = frame["resistant"].astype(int)

    duplicated = int(frame["isolate_id"].duplicated().sum())
    if duplicated:
        raise SystemExit(f"Phenotype file has {duplicated} duplicate isolate rows for this drug. "
                         "Resolve conflicting DST results before training.")
    frame.attrs["dropped_unusable_labels"] = dropped
    return frame[["isolate_id", "resistant"]]


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an isolate x variant feature matrix for amr.py.")
    parser.add_argument("--variants", required=True, type=Path,
                        help="Long-format variant CSV with SAMPLE,CHROM,POS,REF,ALT.")
    parser.add_argument("--loci", type=Path, default=Path("resources/tb_resistance_loci.csv"),
                        help="Gene coordinate table used for annotation and panel selection.")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Directory to create; refuses to overwrite an existing one.")
    parser.add_argument("--panel", choices=["resistance", "genome-wide"], default="resistance",
                        help="'resistance' keeps only variants inside the locus table; "
                             "'genome-wide' keeps every variant passing the prevalence filter.")
    parser.add_argument("--margin", type=int, default=200,
                        help="Bases added either side of each locus, to retain promoter variants.")
    parser.add_argument("--min-prevalence", type=int, default=20,
                        help="Minimum number of isolates carrying a variant for it to be kept. "
                             "Guards against fitting to sequencing singletons.")
    parser.add_argument("--max-prevalence-fraction", type=float, default=0.99,
                        help="Drop variants present in more than this fraction of isolates; "
                             "near-fixed sites carry almost no discriminative signal.")
    parser.add_argument("--metadata", type=Path,
                        help="Output of fetch_sample_metadata.py. Adds real clinical and "
                             "geographic columns and enables the full-vs-baseline comparison.")
    parser.add_argument("--group-by", choices=["isolate", "study"], default="isolate",
                        help="What populates patient_id, the column amr.py keeps intact "
                             "across the train/test split. 'study' uses the ENA study "
                             "accession from --metadata, so isolates from one study "
                             "cannot straddle the split; this is stricter and usually "
                             "lowers measured performance toward a fairer estimate.")
    parser.add_argument("--require-metadata", action="store_true",
                        help="Keep only isolates that have both a country and a specimen "
                             "category, so the two models are compared on identical rows.")
    parser.add_argument("--require-geography", action="store_true",
                        help="Keep only isolates with a known country. A geographic holdout "
                             "needs the holdout column populated for every row, so this is "
                             "the filter to use before geographic_sweep.py.")
    parser.add_argument("--phenotypes", type=Path,
                        help="Optional per-isolate DST table. Without it the matrix is unlabelled.")
    parser.add_argument("--phenotype-sample-column", default="SAMPLE")
    parser.add_argument("--drug-column",
                        help="Column in the phenotype file holding this drug's R/S outcome.")
    parser.add_argument("--antibiotic",
                        help="Antibiotic name written into the output; defaults to --drug-column.")
    parser.add_argument("--organism", default=DEFAULT_ORGANISM)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    if args.phenotypes and not args.drug_column:
        parser.error("--drug-column is required when --phenotypes is supplied.")
    if args.min_prevalence < 1:
        parser.error("--min-prevalence must be at least 1.")
    if not 0.0 < args.max_prevalence_fraction <= 1.0:
        parser.error("--max-prevalence-fraction must be in (0, 1].")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    progress = not args.quiet

    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} already exists. Choose a new output directory so an "
                         "earlier prepared dataset is not silently replaced.")
    if not args.variants.is_file():
        raise SystemExit(f"Variant file not found: {args.variants}")

    loci = LocusIndex.from_csv(args.loci, args.margin) if args.loci.is_file() else None
    if loci is None and args.panel == "resistance":
        raise SystemExit(f"Locus table {args.loci} not found, which the 'resistance' panel requires.")

    if progress:
        print("Pass 1 of 2: counting variant prevalence across isolates.", file=sys.stderr)
    prevalence, gene_of_key, samples = scan_prevalence(
        args.variants, loci, args.panel == "resistance", progress)

    n_samples = len(samples)
    upper_bound = args.max_prevalence_fraction * n_samples
    retained = [key for key, count in prevalence.items()
                if args.min_prevalence <= count <= upper_bound]
    retained.sort(key=lambda key: (gene_of_key.get(key, ""), int(key.split("_", 1)[0])))

    if not retained:
        raise SystemExit("No variants survived the prevalence filter. Lower --min-prevalence.")

    names: dict[str, str] = {}
    used: set[str] = set()
    for key in retained:
        name = column_name(gene_of_key.get(key, ""), key)
        while name in used:          # column names must stay unique after allele shortening
            name = name + "_b"
        used.add(name)
        names[key] = name

    feature_index = {key: position for position, key in enumerate(retained)}
    sample_index = {sample: position for position, sample in enumerate(samples)}

    if progress:
        print(f"Retained {len(retained):,} of {len(prevalence):,} variants across "
              f"{n_samples:,} isolates.", file=sys.stderr)
        print("Pass 2 of 2: filling the isolate x variant matrix.", file=sys.stderr)
    matrix = build_matrix(args.variants, sample_index, feature_index, progress)

    frame = pd.DataFrame(matrix, columns=[names[key] for key in retained])
    frame.insert(0, "isolate_id", samples)
    # No patient identifier exists in this source. Grouping by isolate is the
    # honest fallback; repeated isolates from one patient can then cross the
    # train/test split, which the run report must record as a limitation.
    frame.insert(1, "patient_id", samples)
    frame.insert(2, "organism", args.organism)

    # Real non-genomic covariates, when available. These are the only route to the
    # brief's baseline-versus-full comparison; synthesising them would make the
    # comparison measure nothing.
    clinical_categorical: list[str] = []
    geographic_columns: list[str] = []
    metadata_summary: dict[str, object] = {"metadata_supplied": False}
    if args.metadata:
        if not args.metadata.is_file():
            raise SystemExit(f"Metadata file not found: {args.metadata}")
        meta = pd.read_csv(args.metadata, dtype=str, keep_default_na=False, na_filter=False)
        if "isolate_id" not in meta.columns:
            raise SystemExit("Metadata file has no 'isolate_id' column.")
        wanted = {"country_name": "country", "region": "city_or_region",
                  "specimen_category": "infection_site", "host": "host",
                  "collection_year": "collection_year", "study_accession": "study_accession"}
        available = {src: dest for src, dest in wanted.items() if src in meta.columns}
        meta = meta[["isolate_id", *available]].rename(columns=available)
        meta = meta.drop_duplicates(subset="isolate_id")

        before_join = len(frame)
        frame = frame.merge(meta, on="isolate_id", how="left")
        for column in available.values():
            frame[column] = frame[column].fillna("")

        if args.require_metadata or args.require_geography:
            keep = pd.Series(True, index=frame.index)
            if "country" in frame.columns:
                keep &= frame["country"] != ""
            if args.require_metadata and "infection_site" in frame.columns:
                keep &= frame["infection_site"] != ""
            frame = frame.loc[keep].reset_index(drop=True)

        if args.group_by == "study":
            if "study_accession" not in frame.columns:
                raise SystemExit("--group-by study needs a study_accession column in --metadata.")
            # An isolate with no recorded study keeps its own identity rather than
            # being pooled with every other unknown into one enormous fake group.
            study = frame["study_accession"].where(
                frame["study_accession"] != "", "noStudy_" + frame["isolate_id"])
            frame["patient_id"] = study

        # study_accession is provenance and a grouping key, never a predictor:
        # it identifies the cohort, which is exactly the shortcut to avoid.
        clinical_categorical = [c for c in ("infection_site", "host") if c in frame.columns]
        geographic_columns = [c for c in ("country", "city_or_region") if c in frame.columns]
        metadata_summary = {
            "metadata_supplied": True,
            "grouping": args.group_by,
            "distinct_groups": int(frame["patient_id"].nunique()),
            "metadata_file": str(args.metadata),
            "restricted_to_complete_metadata": bool(args.require_metadata),
            "restricted_to_known_country": bool(args.require_geography),
            "isolates_before_metadata_join": before_join,
            "isolates_after_metadata_join": len(frame),
            "clinical_categorical_columns": clinical_categorical,
            "geographic_columns": geographic_columns,
            "note": ("collection_year is written into the matrix as context but is left out "
                     "of the generated config: it is a temporal covariate, not a clinical "
                     "one, and including it silently would misdescribe the feature blocks."),
        }
        if frame.empty:
            raise SystemExit("No isolate survived the metadata join. Relax --require-metadata.")

    label_summary: dict[str, object] = {"phenotypes_supplied": False}
    antibiotic = args.antibiotic or args.drug_column or ""
    if args.phenotypes:
        phenotypes = load_phenotypes(args.phenotypes, args.phenotype_sample_column, args.drug_column)
        before = len(frame)
        frame = frame.merge(phenotypes, on="isolate_id", how="inner")
        frame.insert(3, "antibiotic", antibiotic)
        counts = frame["resistant"].value_counts().to_dict()
        label_summary = {
            "phenotypes_supplied": True,
            "drug_column": args.drug_column,
            "antibiotic": antibiotic,
            "labels_dropped_as_unusable": phenotypes.attrs.get("dropped_unusable_labels", 0),
            "isolates_before_label_join": before,
            "isolates_with_usable_label": len(frame),
            "resistant": int(counts.get(1, 0)),
            "susceptible": int(counts.get(0, 0)),
        }
        if frame.empty:
            raise SystemExit("No isolate in the variant file matched a usable phenotype. "
                             "Check that the identifier formats agree.")
    else:
        frame.insert(3, "antibiotic", antibiotic)

    args.out_dir.mkdir(parents=True)
    matrix_path = args.out_dir / "isolate_variant_matrix.csv"
    frame.to_csv(matrix_path, index=False)

    dictionary = pd.DataFrame({
        "feature": [names[key] for key in retained],
        "position": [int(key.split("_", 1)[0]) for key in retained],
        "ref": [key.split("_", 2)[1] for key in retained],
        "alt": [key.split("_", 2)[2] for key in retained],
        "gene": [gene_of_key.get(key, "") for key in retained],
        "isolates_carrying": [prevalence[key] for key in retained],
    })
    drug_map = loci.drugs_by_gene() if loci else {}
    note_map = loci.notes_by_gene() if loci else {}
    dictionary["associated_drugs"] = dictionary["gene"].map(drug_map).fillna("")
    dictionary["gene_note"] = dictionary["gene"].map(note_map).fillna("")
    dictionary["carrier_fraction"] = (dictionary["isolates_carrying"] / n_samples).round(5)
    dictionary.to_csv(args.out_dir / "feature_dictionary.csv", index=False)

    genomic_columns = [names[key] for key in retained]
    config = {
        "organism": args.organism,
        "antibiotic": antibiotic,
        "genomic_columns": genomic_columns,
        "organism_column": "organism",
        "antibiotic_column": "antibiotic",
        "target_column": "resistant",
        "sample_column": "isolate_id",
        "group_column": "patient_id",
        "clinical_numeric_columns": [],
        "clinical_categorical_columns": clinical_categorical,
        "geographic_columns": geographic_columns,
        # amr.py only permits the full comparison when genuine clinical AND
        # geographic columns are present, which is exactly the condition here.
        "compare_full": bool(clinical_categorical and geographic_columns),
        "split_method": "group_random",
        "holdout_column": "",
        "holdout_values": [],
        "test_size": 0.2,
        "threshold": 0.5,
        "n_estimators": 300,
        "random_state": 42,
        "permutation_repeats": 0,
    }
    (args.out_dir / "config.generated.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8")

    per_gene = dictionary.groupby("gene")["feature"].count().sort_values(ascending=False)
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "variant_file": str(args.variants),
        "locus_table": str(args.loci) if loci else None,
        "panel": args.panel,
        "locus_margin_bases": args.margin,
        "min_prevalence_isolates": args.min_prevalence,
        "max_prevalence_fraction": args.max_prevalence_fraction,
        "isolates": n_samples,
        "distinct_variants_seen": len(prevalence),
        "features_retained": len(retained),
        "features_per_gene": {str(k): int(v) for k, v in per_gene.items()},
        "metadata": metadata_summary,
        "labels": label_summary,
        "limitations": [
            "No patient identifier exists in the source, so patient_id repeats isolate_id; "
            "repeated isolates from one patient may cross the train/test split.",
            "Variant presence is binary and unfiltered by call quality, depth or allele "
            "fraction, because the source table carries no quality fields.",
            "Locus windows are H37Rv coordinate ranges, not amino-acid level annotations; "
            "a retained variant is not thereby a known resistance-conferring mutation.",
            "No resistance phenotype is derived from the variants themselves.",
        ],
    }
    (args.out_dir / "prep_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {matrix_path} ({len(frame):,} isolates x {len(retained):,} variant features).")
    if not label_summary["phenotypes_supplied"]:
        print("No phenotypes supplied: the matrix has no 'resistant' column and cannot be "
              "used for training. Supply --phenotypes with a documented DST table.")
    else:
        print(f"Labelled isolates: {label_summary['isolates_with_usable_label']:,} "
              f"({label_summary['resistant']:,} resistant / "
              f"{label_summary['susceptible']:,} susceptible).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
