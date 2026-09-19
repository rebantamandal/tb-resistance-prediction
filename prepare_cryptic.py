"""Build a training matrix from the CRyPTIC consortium tables.

This is the CRyPTIC equivalent of `prepare_variants.py`, and it exists because
CRyPTIC closes three gaps the FORUM-TB data leaves open:

  * **Patient identifiers.** Every CRyPTIC sample name carries the subject it
    came from, so repeat samples from one patient can be kept on one side of a
    train/test split. FORUM-TB has no patient field at all.
  * **Genetic lineage.** The project's central explanation is that a model can
    score well by recognising which collection a sample came from, because
    samples from one collection tend to be closely related. Lineage makes that
    testable directly: hold out a whole lineage and see what survives.
  * **One shared laboratory protocol** across all collection sites, so "which
    site" no longer also means "which laboratory method".

It is also a better starting point mechanically. CRyPTIC reports mutations
already annotated with their gene, their protein-level change and a quality
flag, so features come out as `rpoB@S450L` rather than a bare coordinate, and
low-quality calls can be dropped rather than taken at face value.

    python prepare_cryptic.py --drug RIF --out-dir data/cryptic_rif

Run `python prepare_cryptic.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

CHUNK_ROWS = 2_000_000
DEFAULT_ORGANISM = "Mycobacterium tuberculosis"
# site.SITEID.subj.SUBJID.lab.LABID.iso.N
UNIQUEID_PATTERN = re.compile(r"^site\.(?P<site>[^.]+)\.subj\.(?P<subj>.+?)\.lab\.")
VALID_LABELS = {"R": 1, "S": 0}


def parse_uniqueid(series: pd.Series) -> pd.DataFrame:
    """Split CRyPTIC sample names into their site and subject parts."""
    extracted = series.str.extract(UNIQUEID_PATTERN)
    return pd.DataFrame({"site": extracted["site"].fillna(""),
                         "subject": extracted["subj"].fillna("")})


def feature_name(gene: str, mutation: str) -> str:
    """`rpoB@S450L`, kept readable and safe as a column name."""
    gene = re.sub(r"[^A-Za-z0-9_]", "", str(gene))
    mutation = re.sub(r"[^A-Za-z0-9_<>=!*-]", "", str(mutation))
    return f"{gene}@{mutation}"


def load_gene_panel(path: Path) -> set[str]:
    frame = pd.read_csv(path, dtype=str)
    if "gene" not in frame.columns:
        raise SystemExit(f"{path} has no 'gene' column.")
    return {g.strip() for g in frame["gene"] if isinstance(g, str) and g.strip()}


def load_labels(path: Path, drug: str) -> pd.DataFrame:
    """One row per sample with a usable R/S result for this drug.

    Intermediate and unknown outcomes are dropped and counted, never recoded,
    and a sample with conflicting results is discarded rather than arbitrated.
    """
    frame = pd.read_csv(path, compression="infer", dtype=str, low_memory=False)
    for column in ("UNIQUEID", "DRUG", "PHENOTYPE"):
        if column not in frame.columns:
            raise SystemExit(f"{path} has no '{column}' column.")

    subset = frame[frame["DRUG"].str.upper() == drug.upper()].copy()
    if subset.empty:
        available = sorted(frame["DRUG"].dropna().unique())
        raise SystemExit(f"No results for drug '{drug}'. Available: {available}")

    dropped = int((~subset["PHENOTYPE"].isin(VALID_LABELS)).sum())
    subset = subset[subset["PHENOTYPE"].isin(VALID_LABELS)]
    subset["resistant"] = subset["PHENOTYPE"].map(VALID_LABELS)

    grouped = subset.groupby("UNIQUEID")["resistant"].agg(set)
    conflicts = int(sum(1 for values in grouped if len(values) > 1))
    resolved = grouped[grouped.map(len) == 1].map(lambda v: next(iter(v)))

    result = resolved.reset_index()
    result.attrs["dropped_unusable"] = dropped
    result.attrs["conflicts"] = conflicts
    return result


def collect_calls(path: Path, genes: set[str] | None, keep_synonymous: bool,
                  require_pass: bool, wanted: set[str],
                  progress: bool) -> tuple[Counter, dict[str, str], list[tuple[str, str]]]:
    """Read the mutation table once and keep every call that survives filtering.

    A two-pass design would mean decompressing 1.4 GB twice. Restricting to the
    labelled samples and to the resistance-gene panel throws away almost every
    row, so what remains is small enough to hold in memory and reuse: roughly
    fifty mutations per sample rather than the ~1,000 each carries genome-wide.
    """
    prevalence: Counter = Counter()
    gene_of: dict[str, str] = {}
    calls: list[tuple[str, str]] = []
    rows_read = 0

    reader = pd.read_csv(path, compression="infer", chunksize=CHUNK_ROWS, dtype=str,
                         usecols=["UNIQUEID", "GENE", "MUTATION", "IS_SYNONYMOUS",
                                  "IS_FILTER_PASS", "IS_NULL", "IS_HET"],
                         low_memory=False)
    for chunk in reader:
        rows_read += len(chunk)
        chunk = chunk[chunk["UNIQUEID"].isin(wanted)]
        if genes is not None:
            chunk = chunk[chunk["GENE"].isin(genes)]
        if require_pass:
            # Quality flags the FORUM-TB table simply does not have.
            chunk = chunk[chunk["IS_FILTER_PASS"].str.lower() == "true"]
            chunk = chunk[chunk["IS_NULL"].str.lower() != "true"]
            chunk = chunk[chunk["IS_HET"].str.lower() != "true"]
        if not keep_synonymous:
            # A synonymous change does not alter the protein, so it cannot cause
            # resistance; keeping it invites the model to learn lineage markers.
            chunk = chunk[chunk["IS_SYNONYMOUS"].str.lower() != "true"]
        if not chunk.empty:
            names = [feature_name(g, m)
                     for g, m in zip(chunk["GENE"], chunk["MUTATION"])]
            pairs = pd.DataFrame({"s": chunk["UNIQUEID"].to_numpy(), "f": names})
            pairs = pairs.drop_duplicates()
            prevalence.update(pairs["f"].tolist())
            calls.extend(zip(pairs["s"].tolist(), pairs["f"].tolist()))
            for name, gene in zip(names, chunk["GENE"]):
                gene_of.setdefault(name, str(gene))
        if progress:
            print(f"  read {rows_read:,} rows, kept {len(calls):,} calls, "
                  f"{len(prevalence):,} distinct mutations",
                  file=sys.stderr, flush=True)
    return prevalence, gene_of, calls


def build_matrix(calls: list[tuple[str, str]], sample_index: dict[str, int],
                 feature_index: dict[str, int]) -> np.ndarray:
    matrix = np.zeros((len(sample_index), len(feature_index)), dtype=np.uint8)
    for sample, feature in calls:
        column = feature_index.get(feature)
        if column is not None:
            matrix[sample_index[sample], column] = 1
    return matrix


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a matrix from CRyPTIC tables.")
    parser.add_argument("--drug", required=True,
                        help="CRyPTIC three-letter code, e.g. RIF, INH, EMB, PZA, STM.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--mutations", type=Path,
                        default=Path("data/cryptic/MUTATIONS.csv.gz"))
    parser.add_argument("--phenotypes", type=Path,
                        default=Path("data/cryptic/DST_MEASUREMENTS.csv.gz"))
    parser.add_argument("--lineage", type=Path,
                        default=Path("data/cryptic/MYKROBE_LINEAGE.csv.gz"))
    parser.add_argument("--loci", type=Path,
                        default=Path("resources/tb_resistance_loci.csv"))
    parser.add_argument("--panel", choices=["resistance", "genome-wide"],
                        default="resistance")
    parser.add_argument("--min-prevalence", type=int, default=20)
    parser.add_argument("--group-by", choices=["subject", "site", "lineage"],
                        default="subject",
                        help="What fills patient_id, the column amr.py keeps intact "
                             "across the split. 'subject' keeps one patient's samples "
                             "together; 'site' keeps a collection together; 'lineage' "
                             "keeps a genetic family together, which tests whether "
                             "measured accuracy rests on population structure.")
    parser.add_argument("--keep-synonymous", action="store_true",
                        help="Keep silent mutations. Off by default: they cannot cause "
                             "resistance, so they mostly act as lineage markers.")
    parser.add_argument("--no-quality-filter", action="store_true",
                        help="Keep calls that failed quality, are null or heterozygous.")
    parser.add_argument("--antibiotic", help="Name written into the output.")
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    progress = not args.quiet
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} already exists. Choose a new directory.")
    for path in (args.mutations, args.phenotypes):
        if not path.is_file():
            raise SystemExit(f"Not found: {path}")

    if progress:
        print(f"Reading {args.drug} results.", file=sys.stderr)
    labels = load_labels(args.phenotypes, args.drug)
    wanted = set(labels["UNIQUEID"])
    if progress:
        print(f"  {len(wanted):,} samples have a usable result.", file=sys.stderr)

    genes = None
    if args.panel == "resistance":
        if not args.loci.is_file():
            raise SystemExit(f"Locus table not found: {args.loci}")
        genes = load_gene_panel(args.loci)

    if progress:
        print("Reading the mutation table (one pass).", file=sys.stderr)
    prevalence, gene_of, calls = collect_calls(
        args.mutations, genes, args.keep_synonymous,
        not args.no_quality_filter, wanted, progress)

    retained = sorted(f for f, n in prevalence.items() if n >= args.min_prevalence)
    if not retained:
        raise SystemExit("No mutation passed the prevalence filter. Lower it.")
    if progress:
        print(f"Kept {len(retained):,} of {len(prevalence):,} mutations.", file=sys.stderr)

    samples = sorted(wanted)
    sample_index = {s: i for i, s in enumerate(samples)}
    feature_index = {f: i for i, f in enumerate(retained)}
    matrix = build_matrix(calls, sample_index, feature_index)

    frame = pd.DataFrame(matrix, columns=retained)
    frame.insert(0, "isolate_id", samples)
    parts = parse_uniqueid(pd.Series(samples))
    frame.insert(1, "site", parts["site"].to_numpy())
    frame.insert(2, "subject", parts["subject"].to_numpy())

    lineage = pd.Series("", index=range(len(samples)))
    if args.lineage.is_file():
        table = pd.read_csv(args.lineage, compression="infer", dtype=str)
        mapping = dict(zip(table["UNIQUEID"], table["MYKROBE_LINEAGE_NAME_1"].fillna("")))
        lineage = pd.Series([mapping.get(s, "") for s in samples])
    frame.insert(3, "lineage", lineage.to_numpy())

    # An unrecorded grouping value must not pool every unknown into one giant
    # group, so those rows fall back to their own identity.
    source = {"subject": frame["subject"], "site": frame["site"],
              "lineage": frame["lineage"]}[args.group_by]
    group = source.where(source.astype(str).str.strip() != "",
                         f"no{args.group_by}_" + frame["isolate_id"])
    frame.insert(4, "patient_id", group.to_numpy())

    frame.insert(5, "organism", DEFAULT_ORGANISM)
    antibiotic = args.antibiotic or args.drug.upper()
    frame.insert(6, "antibiotic", antibiotic)
    frame = frame.merge(labels.rename(columns={"UNIQUEID": "isolate_id"}),
                        on="isolate_id", how="inner")

    args.out_dir.mkdir(parents=True)
    frame.to_csv(args.out_dir / "isolate_variant_matrix.csv", index=False)

    dictionary = pd.DataFrame({
        "feature": retained,
        "gene": [gene_of.get(f, "") for f in retained],
        "mutation": [f.split("@", 1)[1] if "@" in f else "" for f in retained],
        "samples_carrying": [prevalence[f] for f in retained],
    })
    dictionary.to_csv(args.out_dir / "feature_dictionary.csv", index=False)

    config = {
        "organism": DEFAULT_ORGANISM, "antibiotic": antibiotic,
        "genomic_columns": retained,
        "organism_column": "organism", "antibiotic_column": "antibiotic",
        "target_column": "resistant", "sample_column": "isolate_id",
        "group_column": "patient_id",
        "clinical_numeric_columns": [], "clinical_categorical_columns": [],
        "geographic_columns": [], "compare_full": False,
        "split_method": "group_random", "holdout_column": "", "holdout_values": [],
        "test_size": 0.2, "threshold": 0.5, "n_estimators": 300,
        "random_state": 42, "permutation_repeats": 0,
    }
    (args.out_dir / "config.generated.json").write_text(
        json.dumps(config, indent=2) + "\n", encoding="utf-8")

    counts = frame["resistant"].value_counts().to_dict()
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "CRyPTIC consortium tables",
        "drug": args.drug, "antibiotic": antibiotic,
        "panel": args.panel,
        "quality_filtered": not args.no_quality_filter,
        "synonymous_kept": bool(args.keep_synonymous),
        "min_prevalence": args.min_prevalence,
        "grouped_by": args.group_by,
        "samples": len(frame),
        "resistant": int(counts.get(1, 0)),
        "susceptible": int(counts.get(0, 0)),
        "features": len(retained),
        "distinct_mutations_seen": len(prevalence),
        "distinct_groups": int(frame["patient_id"].nunique()),
        "distinct_sites": int(frame["site"].nunique()),
        "distinct_subjects": int(frame["subject"].nunique()),
        "distinct_lineages": int(frame.loc[frame["lineage"] != "", "lineage"].nunique()),
        "labels_dropped_unusable": labels.attrs.get("dropped_unusable", 0),
        "labels_dropped_conflicting": labels.attrs.get("conflicts", 0),
        "limitations": [
            "Every collection site sits in one country, so site and country remain "
            "linked; what CRyPTIC removes is the differing laboratory protocol, "
            "not the differing patient population.",
            "Lineage is recorded for some sites and not others, so a lineage-grouped "
            "split does not cover the whole cohort.",
        ],
    }
    (args.out_dir / "prep_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {args.out_dir / 'isolate_variant_matrix.csv'} "
          f"({len(frame):,} samples x {len(retained):,} mutations).")
    print(f"  {counts.get(1, 0):,} resistant / {counts.get(0, 0):,} susceptible")
    print(f"  grouped by {args.group_by}: {report['distinct_groups']:,} groups "
          f"({report['distinct_sites']} sites, {report['distinct_subjects']:,} subjects, "
          f"{report['distinct_lineages']} lineages)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
