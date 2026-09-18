"""Recover laboratory drug susceptibility results for ENA/SRA run accessions from BV-BRC.

The variant table carries no resistance phenotype, and a phenotype cannot be
derived from the variants themselves without circularity. BV-BRC (formerly
PATRIC) aggregates per-genome AMR results from the published literature, and
many M. tuberculosis genomes there are named after the sequencing run accession
that produced them, which is what makes the join possible.

Only records whose `evidence` field is "Laboratory Method" are requested. BV-BRC
also stores predictions made by its own classifiers under "Computational
Method"; training a resistance model on those would mean learning another
model's output, so they are excluded at the query and never enter the output.

Run `python fetch_phenotypes.py --help` for options.
"""
from __future__ import annotations

import argparse
import io
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

BVBRC_URL = "https://www.bv-brc.org/api/genome_amr/"
MTB_TAXON = 1773
PAGE_SIZE = 25_000
SELECT_FIELDS = ["genome_id", "genome_name", "antibiotic", "resistant_phenotype",
                 "laboratory_typing_method", "testing_standard", "pubmed", "source"]
ACCESSION_PATTERN = re.compile(r"\b((?:ERR|SRR|DRR)\d{5,})\b")
USER_AGENT = "amr-prediction-research/1.0 (BV-BRC public API)"

# BV-BRC drug names are not consistent with the names used in TB reporting.
DRUG_ALIASES = {
    "rifampin": "RIFAMPICIN",
    "rifampicin": "RIFAMPICIN",
    "isoniazid": "ISONIAZID",
    "ethambutol": "ETHAMBUTOL",
    "pyrazinamide": "PYRAZINAMIDE",
    "streptomycin": "STREPTOMYCIN",
    "amikacin": "AMIKACIN",
    "kanamycin": "KANAMYCIN",
    "capreomycin": "CAPREOMYCIN",
    "ofloxacin": "OFLOXACIN",
    "levofloxacin": "LEVOFLOXACIN",
    "moxifloxacin": "MOXIFLOXACIN",
    "ciprofloxacin": "CIPROFLOXACIN",
    "ethionamide": "ETHIONAMIDE",
    "prothionamide": "PROTHIONAMIDE",
    "cycloserine": "CYCLOSERINE",
    "para-aminosalicylic acid": "PARA_AMINOSALICYLIC_ACID",
    "bedaquiline": "BEDAQUILINE",
    "clofazimine": "CLOFAZIMINE",
    "linezolid": "LINEZOLID",
    "delamanid": "DELAMANID",
    "rifabutin": "RIFABUTIN",
}

# Only unambiguous outcomes become labels. Intermediate is a real laboratory
# result, but it is neither of the two classes being modelled, so it is dropped
# and counted rather than folded into either one.
PHENOTYPE_MAP = {
    "resistant": "R",
    "susceptible": "S",
    "sensitive": "S",
}
DROP_PHENOTYPES = {"intermediate", "not defined", "inconclusive", "indeterminate", ""}


def request_page(offset: int, timeout: int, retries: int, pause: float) -> pd.DataFrame:
    """One page of laboratory-evidence AMR records, with backoff."""
    query = (f"and(eq(taxon_id,{MTB_TAXON}),eq(evidence,%22Laboratory%20Method%22))"
             f"&select({','.join(SELECT_FIELDS)})"
             f"&limit({PAGE_SIZE},{offset})")
    url = f"{BVBRC_URL}?{query}"

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            url, headers={"User-Agent": USER_AGENT, "Accept": "text/tsv"})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
            if not text.strip():
                return pd.DataFrame(columns=SELECT_FIELDS)
            frame = pd.read_csv(io.StringIO(text), sep="\t", dtype=str,
                                keep_default_na=False, na_filter=False)
            # BV-BRC quotes its TSV values; strip the quoting rather than trusting it.
            for column in frame.columns:
                frame[column] = frame[column].str.strip().str.strip('"')
            return frame
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as error:
            last_error = error
            if attempt < retries:
                delay = pause * (2 ** (attempt - 1))
                print(f"    request failed ({error}); retrying in {delay:.1f}s",
                      file=sys.stderr, flush=True)
                time.sleep(delay)
    raise SystemExit(f"BV-BRC request failed after {retries} attempts: {last_error}")


def fetch_genome_index(cache_dir: Path, timeout: int, retries: int, pause: float,
                       progress: bool) -> dict[str, str]:
    """Map BV-BRC genome_id to its registered SRA/ENA accession.

    Most M. tuberculosis genomes in BV-BRC are named after a laboratory strain
    rather than a run accession, so parsing the name alone finds only a fraction
    of them. The genome table records the originating accession explicitly,
    which recovers far more of the cohort.
    """
    index: dict[str, str] = {}
    for page in range(20):
        offset = page * PAGE_SIZE
        cached = cache_dir / f"genome_{page:03d}.tsv"
        if cached.is_file():
            frame = pd.read_csv(cached, sep="\t", dtype=str,
                                keep_default_na=False, na_filter=False)
        else:
            query = (f"and(eq(taxon_id,{MTB_TAXON}),ne(sra_accession,%22%22))"
                     f"&select(genome_id,sra_accession)&limit({PAGE_SIZE},{offset})")
            url = f"https://www.bv-brc.org/api/genome/?{query}"
            last_error: Exception | None = None
            frame = None
            for attempt in range(1, retries + 1):
                request = urllib.request.Request(
                    url, headers={"User-Agent": USER_AGENT, "Accept": "text/tsv"})
                try:
                    with urllib.request.urlopen(request, timeout=timeout) as response:
                        text = response.read().decode("utf-8", errors="replace")
                    frame = (pd.read_csv(io.StringIO(text), sep="\t", dtype=str,
                                         keep_default_na=False, na_filter=False)
                             if text.strip() else pd.DataFrame(columns=["genome_id",
                                                                        "sra_accession"]))
                    for column in frame.columns:
                        frame[column] = frame[column].str.strip().str.strip('"')
                    break
                except (urllib.error.URLError, urllib.error.HTTPError,
                        TimeoutError, OSError) as error:
                    last_error = error
                    if attempt < retries:
                        time.sleep(pause * (2 ** (attempt - 1)))
            if frame is None:
                raise SystemExit(f"BV-BRC genome request failed: {last_error}")
            frame.to_csv(cached, sep="\t", index=False)
            time.sleep(pause)

        if frame.empty:
            break
        for genome_id, accession in zip(frame["genome_id"], frame["sra_accession"]):
            if genome_id and accession:
                index[genome_id] = accession
        if progress:
            print(f"  genome page {page}: {len(frame)} rows "
                  f"({len(index):,} accessions)", file=sys.stderr, flush=True)
        if len(frame) < PAGE_SIZE:
            break
    return index


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch laboratory DST results for run accessions from BV-BRC.")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output phenotype CSV; refuses to overwrite.")
    parser.add_argument("--accessions", type=Path,
                        help="Optional matrix/list to report overlap against and restrict to.")
    parser.add_argument("--id-column", default="isolate_id")
    parser.add_argument("--cache", type=Path,
                        help="Directory for raw pages, making the run resumable.")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--pause", type=float, default=1.0)
    parser.add_argument("--max-pages", type=int, default=40)
    parser.add_argument("--quiet", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    progress = not args.quiet

    if args.out.exists():
        raise SystemExit(f"{args.out} already exists. Choose a new filename.")

    cache_dir = args.cache or args.out.parent / ".phenotype_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    pages: list[pd.DataFrame] = []
    for page in range(args.max_pages):
        offset = page * PAGE_SIZE
        cached = cache_dir / f"page_{page:03d}.tsv"
        if cached.is_file():
            frame = pd.read_csv(cached, sep="\t", dtype=str,
                                keep_default_na=False, na_filter=False)
            if progress:
                print(f"  page {page}: cached ({len(frame)} rows)", file=sys.stderr, flush=True)
        else:
            frame = request_page(offset, args.timeout, args.retries, args.pause)
            frame.to_csv(cached, sep="\t", index=False)
            if progress:
                print(f"  page {page}: {len(frame)} rows", file=sys.stderr, flush=True)
            time.sleep(args.pause)
        if frame.empty:
            break
        pages.append(frame)
        if len(frame) < PAGE_SIZE:
            break

    if not pages:
        raise SystemExit("BV-BRC returned no laboratory-evidence records.")
    records = pd.concat(pages, ignore_index=True)
    total_records = len(records)

    # Two routes to a run accession: the genome table's registered SRA accession,
    # and the accession embedded in the genome name. The registered accession is
    # authoritative, so it is tried first and the name is only a fallback.
    if progress:
        print("Fetching genome accession index.", file=sys.stderr)
    genome_index = fetch_genome_index(cache_dir, args.timeout, args.retries,
                                      args.pause, progress)
    registered = records["genome_id"].map(genome_index).fillna("")
    registered = registered.where(registered.str.match(ACCESSION_PATTERN.pattern.strip("\\b")),
                                  "")
    from_name = records["genome_name"].str.extract(ACCESSION_PATTERN, expand=False).fillna("")
    records["isolate_id"] = registered.where(registered != "", from_name)

    matched_via_genome_table = int((registered != "").sum())
    matched_via_name_only = int(((registered == "") & (from_name != "")).sum())
    records = records[records["isolate_id"] != ""]
    with_accession = len(records)

    records["drug"] = records["antibiotic"].str.strip().str.lower().map(DRUG_ALIASES)
    unmapped = sorted(set(records.loc[records["drug"].isna(), "antibiotic"].str.lower()))
    records = records.dropna(subset=["drug"])

    raw_phenotype = records["resistant_phenotype"].str.strip().str.lower()
    records["label"] = raw_phenotype.map(PHENOTYPE_MAP)
    dropped_ambiguous = int(records["label"].isna().sum())
    dropped_values = Counter(raw_phenotype[records["label"].isna()])
    records = records.dropna(subset=["label"])

    if args.accessions is not None:
        if args.accessions.suffix.lower() == ".csv":
            wanted = pd.read_csv(args.accessions, usecols=[args.id_column], dtype=str)
            wanted = set(wanted[args.id_column].dropna())
        else:
            wanted = {line.strip() for line in
                      args.accessions.read_text(encoding="utf-8").splitlines() if line.strip()}
    else:
        wanted = None

    # One isolate and drug can appear more than once, from different publications.
    # Agreeing duplicates collapse; genuinely conflicting results are discarded,
    # because silently preferring one laboratory over another would be arbitrary.
    grouped = records.groupby(["isolate_id", "drug"])["label"].agg(set)
    conflicts = int(sum(1 for values in grouped if len(values) > 1))
    resolved = grouped[grouped.map(len) == 1].map(lambda values: next(iter(values)))
    wide = resolved.unstack("drug")
    wide.index.name = "SAMPLE"
    wide = wide.reset_index().fillna("")

    overlap_summary: dict[str, object] = {}
    if wanted is not None:
        overlapping = wide[wide["SAMPLE"].isin(wanted)].reset_index(drop=True)
        overlap_summary = {
            "accessions_in_variant_file": len(wanted),
            "accessions_with_any_laboratory_dst": len(overlapping),
        }
        per_drug = {}
        for drug in sorted(c for c in overlapping.columns if c != "SAMPLE"):
            counts = Counter(v for v in overlapping[drug] if v)
            if counts:
                per_drug[drug] = {"labelled": sum(counts.values()),
                                  "resistant": counts.get("R", 0),
                                  "susceptible": counts.get("S", 0)}
        overlap_summary["per_drug_within_variant_file"] = dict(
            sorted(per_drug.items(), key=lambda kv: -kv[1]["labelled"]))
        wide = overlapping

    args.out.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(args.out, index=False)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "BV-BRC genome_amr API, taxon_id=1773, evidence='Laboratory Method' only",
        "computational_predictions_excluded": True,
        "records_retrieved": total_records,
        "records_with_run_accession": with_accession,
        "matched_via_genome_table_sra_accession": matched_via_genome_table,
        "matched_via_genome_name_only": matched_via_name_only,
        "records_dropped_unmappable_drug_name": unmapped,
        "records_dropped_ambiguous_phenotype": dropped_ambiguous,
        "dropped_phenotype_values": dict(dropped_values.most_common()),
        "isolate_drug_pairs_with_conflicting_results": conflicts,
        "overlap": overlap_summary,
        "limitations": [
            "BV-BRC aggregates results from many publications using different DST "
            "methods, critical concentrations and standards; they are not harmonised here.",
            "Only genomes whose BV-BRC name embeds the run accession can be joined, so "
            "coverage is a floor, not the true extent of published DST for this cohort.",
            "Conflicting results for the same isolate and drug are discarded rather than "
            "arbitrated; review them manually if coverage matters more than caution.",
            "Pyrazinamide DST in particular is known to be poorly reproducible.",
        ],
    }
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {args.out} ({len(wide):,} isolates with at least one laboratory DST result).")
    if overlap_summary.get("per_drug_within_variant_file"):
        for drug, counts in overlap_summary["per_drug_within_variant_file"].items():
            print(f"  {drug:26s} {counts['labelled']:5d} labelled "
                  f"({counts['resistant']} R / {counts['susceptible']} S)")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
