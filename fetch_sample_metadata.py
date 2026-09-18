"""Retrieve real per-isolate metadata for ENA/SRA run accessions from the ENA portal API.

The variant table identifies isolates by sequencing run accession (ERR.../SRR...).
Those runs have publicly registered sample metadata: the country the isolate was
collected in, the specimen it was isolated from, the collection date, and the
originating study. This script fetches that metadata and nothing else. It does
not synthesise, impute or infer any field: an isolate whose submitter never
registered a country simply has no country here, and is reported as such.

The ENA portal API is public and needs no account or key.

Fields retrieved and how they map to the project's feature categories:

  country, region        -> geographic features
  isolation_source       -> clinical feature (the specimen / infection site)
  host                   -> clinical context
  collection_year        -> temporal context, useful for a time-based holdout
  study_accession        -> provenance; also a candidate grouping variable

Run `python fetch_sample_metadata.py --help` for options.
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ENA_SEARCH_URL = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FIELDS = [
    "run_accession",
    "sample_accession",
    "country",
    "location",
    "collection_date",
    "isolation_source",
    "host",
    "scientific_name",
    "study_accession",
    "center_name",
]
USER_AGENT = "amr-prediction-research/1.0 (ENA portal API; public metadata only)"


def read_accessions(path: Path, column: str) -> list[str]:
    """Read the isolate identifiers, from either a prepared matrix or a plain list."""
    if path.suffix.lower() in {".csv", ".tsv"}:
        separator = "\t" if path.suffix.lower() == ".tsv" else ","
        # Only the identifier column is needed; the matrix itself can be large.
        frame = pd.read_csv(path, sep=separator, usecols=lambda c: c == column,
                            dtype=str, keep_default_na=False, na_filter=False)
        if column not in frame.columns:
            raise SystemExit(f"{path} has no column '{column}'.")
        values = frame[column].tolist()
    else:
        values = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]

    accessions = sorted({v.strip() for v in values if v.strip()})
    if not accessions:
        raise SystemExit(f"No accessions found in {path}.")
    return accessions


def fetch_batch(accessions: list[str], timeout: int, retries: int,
                pause: float) -> pd.DataFrame:
    """One POST to the ENA portal API, with backoff on transient failures."""
    payload = urllib.parse.urlencode({
        "result": "read_run",
        "includeAccessions": ",".join(accessions),
        "fields": ",".join(ENA_FIELDS),
        "format": "tsv",
        "limit": "0",
    }).encode("utf-8")

    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        request = urllib.request.Request(
            ENA_SEARCH_URL, data=payload,
            headers={"User-Agent": USER_AGENT,
                     "Content-Type": "application/x-www-form-urlencoded"})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                text = response.read().decode("utf-8", errors="replace")
            if not text.strip():
                return pd.DataFrame(columns=ENA_FIELDS)
            return pd.read_csv(io.StringIO(text), sep="\t", dtype=str,
                               keep_default_na=False, na_filter=False)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as error:
            last_error = error
            if attempt < retries:
                delay = pause * (2 ** (attempt - 1))
                print(f"    request failed ({error}); retrying in {delay:.1f}s",
                      file=sys.stderr, flush=True)
                time.sleep(delay)
    raise SystemExit(f"ENA request failed after {retries} attempts: {last_error}")


class Normalizer:
    """Applies the editable rules in resources/metadata_normalization.json.

    The rules only recognise explicit null placeholders and group verbatim
    specimen strings. No rule invents a value a submitter did not record: ENA
    writes the literal string 'missing' into these fields, and treating that as
    data would manufacture a category that means nothing.
    """

    def __init__(self, rules: dict) -> None:
        self.null_tokens = {str(t).strip().lower() for t in rules.get("null_tokens", [])}
        self.facility_values = {str(t).strip().lower() for t in rules.get("facility_values", [])}
        self.specimen_categories = {str(k).strip().lower(): str(v)
                                    for k, v in rules.get("specimen_categories", {}).items()}

    @classmethod
    def load(cls, path: Path | None) -> "Normalizer":
        if path is None or not path.is_file():
            return cls({})
        return cls(json.loads(path.read_text(encoding="utf-8")))

    def is_null(self, value: str) -> bool:
        return (value or "").strip().lower() in self.null_tokens

    def is_facility(self, value: str) -> bool:
        return (value or "").strip().lower() in self.facility_values

    def specimen_category(self, value: str) -> str:
        cleaned = (value or "").strip().lower()
        if not cleaned or cleaned in self.null_tokens:
            return ""
        return self.specimen_categories.get(cleaned, "other_recorded")


def split_country(value: str) -> tuple[str, str]:
    """ENA records country as 'Country: region, detail'. Split without inventing detail."""
    value = (value or "").strip()
    if not value:
        return "", ""
    if ":" in value:
        country, region = value.split(":", 1)
        return country.strip(), region.strip()
    return value, ""


def collection_year(value: str) -> str:
    """Take the year from a partial ISO date, leaving anything unparseable blank."""
    value = (value or "").strip()
    if len(value) >= 4 and value[:4].isdigit():
        year = int(value[:4])
        if 1900 <= year <= datetime.now(timezone.utc).year:
            return str(year)
    return ""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fetch public ENA metadata for sequencing run accessions.")
    parser.add_argument("--accessions", required=True, type=Path,
                        help="Prepared matrix CSV, or a text file with one accession per line.")
    parser.add_argument("--id-column", default="isolate_id",
                        help="Identifier column, when --accessions is a CSV.")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output CSV; refuses to overwrite an existing file.")
    parser.add_argument("--cache", type=Path,
                        help="Directory for per-batch responses, making the run resumable. "
                             "Defaults to a .cache folder beside --out.")
    parser.add_argument("--batch-size", type=int, default=400,
                        help="Accessions per request.")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=4)
    parser.add_argument("--pause", type=float, default=1.0,
                        help="Seconds between requests, to stay courteous to a public service.")
    parser.add_argument("--normalization", type=Path,
                        default=Path("resources/metadata_normalization.json"),
                        help="Editable rules for null placeholders and specimen grouping.")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if args.batch_size < 1 or args.batch_size > 1000:
        parser.error("--batch-size must be between 1 and 1000.")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    progress = not args.quiet

    if args.out.exists():
        raise SystemExit(f"{args.out} already exists. Choose a new filename so an earlier "
                         "metadata pull is not silently replaced.")
    if not args.accessions.is_file():
        raise SystemExit(f"Accession source not found: {args.accessions}")

    accessions = read_accessions(args.accessions, args.id_column)
    cache_dir = args.cache or args.out.parent / ".metadata_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    batches = [accessions[i:i + args.batch_size]
               for i in range(0, len(accessions), args.batch_size)]
    if progress:
        print(f"Fetching metadata for {len(accessions):,} accessions "
              f"in {len(batches)} batches.", file=sys.stderr)

    frames: list[pd.DataFrame] = []
    for number, batch in enumerate(batches, start=1):
        cached = cache_dir / f"batch_{number:04d}.tsv"
        if cached.is_file():
            frames.append(pd.read_csv(cached, sep="\t", dtype=str,
                                      keep_default_na=False, na_filter=False))
            if progress:
                print(f"  batch {number}/{len(batches)}: cached", file=sys.stderr, flush=True)
            continue

        frame = fetch_batch(batch, args.timeout, args.retries, args.pause)
        frame.to_csv(cached, sep="\t", index=False)
        frames.append(frame)
        if progress:
            print(f"  batch {number}/{len(batches)}: {len(frame)} records",
                  file=sys.stderr, flush=True)
        if number < len(batches):
            time.sleep(args.pause)

    fetched = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=ENA_FIELDS)
    for field in ENA_FIELDS:
        if field not in fetched.columns:
            fetched[field] = ""
    fetched = fetched.drop_duplicates(subset="run_accession")

    country_parts = fetched["country"].map(split_country)
    fetched["country_raw"] = [c for c, _ in country_parts]
    fetched["region"] = [r for _, r in country_parts]
    fetched["collection_year"] = fetched["collection_date"].map(collection_year)

    normalizer = Normalizer.load(args.normalization)
    # A facility name recorded in the country field is preserved verbatim in
    # collection_site rather than being resolved to a country by guesswork.
    fetched["collection_site"] = [v if normalizer.is_facility(v) else ""
                                  for v in fetched["country_raw"]]
    fetched["country_name"] = ["" if (normalizer.is_null(v) or normalizer.is_facility(v)) else v
                               for v in fetched["country_raw"]]
    fetched["region"] = ["" if normalizer.is_null(v) else v for v in fetched["region"]]
    fetched["specimen_category"] = fetched["isolation_source"].map(normalizer.specimen_category)
    fetched["host"] = ["" if normalizer.is_null(v) else v for v in fetched["host"]]

    # Left join onto the requested accessions so isolates ENA has no record for
    # stay visible as blanks rather than silently disappearing.
    result = pd.DataFrame({"isolate_id": accessions}).merge(
        fetched, left_on="isolate_id", right_on="run_accession", how="left")
    result = result.fillna("")

    columns = ["isolate_id", "sample_accession", "country_name", "region",
               "specimen_category", "collection_year", "country_raw", "collection_site",
               "collection_date", "isolation_source", "host",
               "scientific_name", "study_accession", "center_name"]
    result[columns].to_csv(args.out, index=False)

    matched = int((result["run_accession"] != "").sum())
    with_country = int((result["country_name"] != "").sum())
    with_source = int((result["specimen_category"] != "").sum())
    with_year = int((result["collection_year"] != "").sum())

    country_counts = Counter(c for c in result["country_name"] if c)
    specimen_counts = Counter(s for s in result["specimen_category"] if s)
    facility_count = int((result["collection_site"] != "").sum())
    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "ENA portal API (https://www.ebi.ac.uk/ena/portal/api), result=read_run",
        "accessions_requested": len(accessions),
        "accessions_matched_in_ena": matched,
        "normalization_rules": str(args.normalization) if args.normalization else None,
        "coverage_after_normalization": {
            "country": with_country,
            "specimen_category": with_source,
            "collection_year": with_year,
            "facility_name_in_country_field": facility_count,
        },
        "distinct_countries": len(country_counts),
        "isolates_per_country": dict(country_counts.most_common()),
        "isolates_per_specimen_category": dict(specimen_counts.most_common()),
        "limitations": [
            "Coverage is whatever the original submitter registered; SRA-sourced runs "
            "are frequently missing country and isolation_source entirely.",
            "Country is the collection country recorded for the sample, not the patient's "
            "residence or country of infection.",
            "ENA free-text fields are not harmonised; isolation_source in particular needs "
            "manual grouping into a small set of categories before use as a feature.",
            "Missing values are left blank and are never imputed here.",
        ],
    }
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Wrote {args.out} ({matched:,} of {len(accessions):,} accessions matched in ENA).")
    print(f"  usable country: {with_country:,}  specimen category: {with_source:,}  "
          f"collection year: {with_year:,}")
    print(f"  {len(country_counts)} distinct countries; top: "
          f"{', '.join(f'{c} ({n})' for c, n in country_counts.most_common(5))}")
    if facility_count:
        print(f"  {facility_count} isolates record a facility name instead of a country; "
              "left blank and preserved in collection_site.")
    print(f"Coverage report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
