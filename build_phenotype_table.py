"""Combine laboratory DST results from several sources into one phenotype table.

Three sources are supported:

  BV-BRC   laboratory-evidence AMR records, joined by the run accession embedded
           in the genome name (see fetch_phenotypes.py).
  CRyPTIC  the consortium's reuse table, which reports binary phenotypes derived
           from a harmonised broth microdilution plate, along with a quality
           grade per measurement.
  NCBI     Pathogen Detection metadata, whose AST_phenotypes field carries the
           submitter-reported susceptibility results for each sequencing run.

Where sources overlap, CRyPTIC is preferred, then NCBI, then BV-BRC: CRyPTIC's
measurements come from one standardised assay with an explicit quality grade,
whereas the other two aggregate many laboratories and standards. Disagreements
are counted and reported for every overlapping pair rather than hidden, because
a high disagreement rate bounds how well any model can possibly score.

Every output label keeps a provenance column naming the source it came from.

Run `python build_phenotype_table.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

# CRyPTIC uses three-letter drug codes in its column names.
CRYPTIC_DRUGS = {
    "AMI": "AMIKACIN",
    "BDQ": "BEDAQUILINE",
    "CFZ": "CLOFAZIMINE",
    "DLM": "DELAMANID",
    "EMB": "ETHAMBUTOL",
    "ETH": "ETHIONAMIDE",
    "INH": "ISONIAZID",
    "KAN": "KANAMYCIN",
    "LEV": "LEVOFLOXACIN",
    "LZD": "LINEZOLID",
    "MXF": "MOXIFLOXACIN",
    "RIF": "RIFAMPICIN",
    "RFB": "RIFABUTIN",
}
VALID_LABELS = {"R", "S"}


def load_cryptic(path: Path, min_quality: set[str]) -> pd.DataFrame:
    """Long-format (isolate, drug, label) from the CRyPTIC reuse table."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
    if "ENA_RUN" not in frame.columns:
        raise SystemExit(f"{path} has no ENA_RUN column; is this the CRyPTIC reuse table?")

    rows = []
    for code, drug in CRYPTIC_DRUGS.items():
        label_column = f"{code}_BINARY_PHENOTYPE"
        quality_column = f"{code}_PHENOTYPE_QUALITY"
        if label_column not in frame.columns:
            continue
        subset = pd.DataFrame({
            "isolate_id": frame["ENA_RUN"],
            "drug": drug,
            "label": frame[label_column].str.strip().str.upper(),
            "quality": (frame[quality_column].str.strip().str.upper()
                        if quality_column in frame.columns else ""),
        })
        subset = subset[subset["label"].isin(VALID_LABELS)]
        if min_quality and quality_column in frame.columns:
            subset = subset[subset["quality"].isin(min_quality)]
        rows.append(subset[["isolate_id", "drug", "label"]])

    if not rows:
        return pd.DataFrame(columns=["isolate_id", "drug", "label"])
    combined = pd.concat(rows, ignore_index=True)
    combined["source"] = "CRyPTIC"
    return combined


NCBI_DRUG_ALIASES = {
    "rifampin": "RIFAMPICIN", "rifampicin": "RIFAMPICIN", "isoniazid": "ISONIAZID",
    "ethambutol": "ETHAMBUTOL", "pyrazinamide": "PYRAZINAMIDE",
    "streptomycin": "STREPTOMYCIN", "amikacin": "AMIKACIN", "kanamycin": "KANAMYCIN",
    "capreomycin": "CAPREOMYCIN", "ofloxacin": "OFLOXACIN",
    "levofloxacin": "LEVOFLOXACIN", "moxifloxacin": "MOXIFLOXACIN",
    "ciprofloxacin": "CIPROFLOXACIN", "ethionamide": "ETHIONAMIDE",
    "cycloserine": "CYCLOSERINE", "bedaquiline": "BEDAQUILINE",
    "clofazimine": "CLOFAZIMINE", "linezolid": "LINEZOLID",
    "clarithromycin": "CLARITHROMYCIN", "rifabutin": "RIFABUTIN",
}


def load_ncbi(path: Path) -> pd.DataFrame:
    """Long-format (isolate, drug, label) from the NCBI Pathogen Detection metadata.

    AST_phenotypes is a comma-separated 'drug=value' list, and one row can list
    several sequencing runs, so the phenotype is attached to each run named.
    Intermediate ('I') is dropped rather than folded into either class.
    """
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                        na_filter=False, low_memory=False)
    for needed in ("Run", "AST_phenotypes"):
        if needed not in frame.columns:
            raise SystemExit(f"{path} has no {needed} column; is this the "
                             "NCBI Pathogen Detection metadata TSV?")

    rows = []
    for runs, phenotypes in zip(frame["Run"], frame["AST_phenotypes"]):
        phenotypes = str(phenotypes).strip()
        if not phenotypes or phenotypes.upper() == "NULL":
            continue
        accessions = [r.strip() for r in str(runs).split(",") if r.strip()]
        if not accessions:
            continue
        for token in phenotypes.split(","):
            if "=" not in token:
                continue
            drug_raw, value = token.split("=", 1)
            drug = NCBI_DRUG_ALIASES.get(drug_raw.strip().lower())
            label = value.strip().upper()
            if drug is None or label not in VALID_LABELS:
                continue
            for accession in accessions:
                rows.append((accession, drug, label))

    result = pd.DataFrame(rows, columns=["isolate_id", "drug", "label"])
    result = result.drop_duplicates()
    result["source"] = "NCBI"
    return result


def load_bvbrc(path: Path) -> pd.DataFrame:
    """Long-format (isolate, drug, label) from the wide BV-BRC output."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False, na_filter=False)
    if "SAMPLE" not in frame.columns:
        raise SystemExit(f"{path} has no SAMPLE column; is this fetch_phenotypes.py output?")
    drugs = [c for c in frame.columns if c != "SAMPLE"]
    melted = frame.melt(id_vars="SAMPLE", value_vars=drugs,
                        var_name="drug", value_name="label")
    melted = melted.rename(columns={"SAMPLE": "isolate_id"})
    melted["label"] = melted["label"].str.strip().str.upper()
    melted = melted[melted["label"].isin(VALID_LABELS)]
    melted["source"] = "BV-BRC"
    return melted[["isolate_id", "drug", "label", "source"]]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge laboratory DST results into one phenotype table.")
    parser.add_argument("--bvbrc", type=Path, help="Output of fetch_phenotypes.py.")
    parser.add_argument("--ncbi", type=Path,
                        help="NCBI Pathogen Detection metadata TSV (has AST_phenotypes).")
    parser.add_argument("--cryptic", type=Path, help="CRyPTIC reuse table CSV.")
    parser.add_argument("--cryptic-quality", default="HIGH,MEDIUM",
                        help="Comma-separated CRyPTIC quality grades to accept, "
                             "or 'ANY' to accept all.")
    parser.add_argument("--accessions", type=Path,
                        help="Restrict output to the isolates in this matrix or list.")
    parser.add_argument("--id-column", default="isolate_id")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    if not args.bvbrc and not args.cryptic and not args.ncbi:
        parser.error("Supply at least one of --bvbrc, --cryptic or --ncbi.")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.out.exists():
        raise SystemExit(f"{args.out} already exists. Choose a new filename.")

    quality = set() if args.cryptic_quality.strip().upper() == "ANY" else {
        q.strip().upper() for q in args.cryptic_quality.split(",") if q.strip()}

    sources: list[pd.DataFrame] = []
    source_counts: dict[str, int] = {}
    if args.cryptic:
        if not args.cryptic.is_file():
            raise SystemExit(f"CRyPTIC table not found: {args.cryptic}")
        cryptic = load_cryptic(args.cryptic, quality)
        source_counts["CRyPTIC"] = len(cryptic)
        sources.append(cryptic)
    if args.ncbi:
        if not args.ncbi.is_file():
            raise SystemExit(f"NCBI metadata not found: {args.ncbi}")
        ncbi = load_ncbi(args.ncbi)
        source_counts["NCBI"] = len(ncbi)
        sources.append(ncbi)
    if args.bvbrc:
        if not args.bvbrc.is_file():
            raise SystemExit(f"BV-BRC table not found: {args.bvbrc}")
        bvbrc = load_bvbrc(args.bvbrc)
        source_counts["BV-BRC"] = len(bvbrc)
        sources.append(bvbrc)

    combined = pd.concat(sources, ignore_index=True)

    if args.accessions is not None:
        if args.accessions.suffix.lower() == ".csv":
            wanted = set(pd.read_csv(args.accessions, usecols=[args.id_column],
                                     dtype=str)[args.id_column].dropna())
        else:
            wanted = {line.strip() for line in
                      args.accessions.read_text(encoding="utf-8").splitlines() if line.strip()}
        combined = combined[combined["isolate_id"].isin(wanted)]

    # Measure agreement before resolving, so label quality is visible.
    pivot = combined.pivot_table(index=["isolate_id", "drug"], columns="source",
                                 values="label", aggfunc="first")
    # Agreement is measured for every pair of sources that actually overlaps.
    # A low rate is itself a finding: it bounds how well any model can score,
    # because the labels themselves would then be partly wrong.
    agreement: dict[str, object] = {}
    present = [c for c in pivot.columns]
    for i, left in enumerate(present):
        for right in present[i + 1:]:
            both = pivot.dropna(subset=[left, right])
            if both.empty:
                agreement[f"{left} vs {right}"] = {"covered_by_both": 0}
                continue
            agree = int((both[left] == both[right]).sum())
            agreement[f"{left} vs {right}"] = {
                "covered_by_both": len(both),
                "agreeing": agree,
                "disagreeing": len(both) - agree,
                "agreement_rate": round(agree / len(both), 4),
            }

    # CRyPTIC wins ties: one harmonised assay with a quality grade beats a
    # mixture of laboratories and standards.
    priority = {"CRyPTIC": 0, "NCBI": 1, "BV-BRC": 2}
    combined["_rank"] = combined["source"].map(priority).fillna(9)
    combined = combined.sort_values("_rank").drop_duplicates(
        subset=["isolate_id", "drug"], keep="first")

    labels = combined.pivot(index="isolate_id", columns="drug", values="label")
    provenance = combined.pivot(index="isolate_id", columns="drug", values="source")
    provenance.columns = [f"{c}_SOURCE" for c in provenance.columns]

    result = labels.join(provenance).fillna("")
    result.index.name = "SAMPLE"
    result = result.reset_index()

    drug_columns = sorted(c for c in labels.columns)
    ordered = ["SAMPLE"] + drug_columns + [f"{d}_SOURCE" for d in drug_columns]
    result = result[ordered]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)

    per_drug = {}
    for drug in drug_columns:
        counts = Counter(v for v in result[drug] if v)
        if counts:
            per_drug[drug] = {"labelled": sum(counts.values()),
                              "resistant": counts.get("R", 0),
                              "susceptible": counts.get("S", 0),
                              "resistant_fraction": round(
                                  counts.get("R", 0) / sum(counts.values()), 4)}

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources_used": source_counts,
        "cryptic_quality_grades_accepted": sorted(quality) if quality else "ANY",
        "tie_break": "Where sources overlap, CRyPTIC is preferred, then NCBI, then BV-BRC",
        "cross_source_agreement": agreement,
        "isolates_with_any_label": len(result),
        "per_drug": dict(sorted(per_drug.items(), key=lambda kv: -kv[1]["labelled"])),
        "limitations": [
            "Labels come from different laboratories, DST methods and critical "
            "concentrations; they are not harmonised beyond the binary R/S call.",
            "A per-drug source column is included so any run can be repeated on a "
            "single source if mixing them proves unsatisfactory.",
            "Pyrazinamide and ethambutol DST are known to be poorly reproducible; "
            "treat their labels, and any model trained on them, with more caution.",
        ],
    }
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {args.out} ({len(result):,} isolates).")
    for drug, counts in report["per_drug"].items():
        print(f"  {drug:26s} {counts['labelled']:5d} labelled  "
              f"{counts['resistant']:5d} R / {counts['susceptible']:5d} S  "
              f"({counts['resistant_fraction']:.1%} resistant)")
    for pair, stats in agreement.items():
        if stats.get("covered_by_both"):
            print(f"Agreement {pair}: {stats['agreeing']}/{stats['covered_by_both']} "
                  f"({stats['agreement_rate']:.1%})")
        else:
            print(f"Agreement {pair}: no shared isolate and drug, cannot be compared.")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
