"""Hold out each country in turn and record how far performance falls.

A random train/test split answers "can the model predict resistance in an
isolate whose relatives it has already seen?". Deployment asks something harder:
"can it predict for a patient in a setting that contributed nothing to
training?". Those are different questions and, on pooled public TB data, they
get very different answers.

This script runs one `amr.py` geographic holdout per country and collects the
results into a single table, so the gap between the two questions is visible at
a glance rather than buried across run directories.

Countries whose held-out set contains only one class are reported as undefined
rather than being silently dropped: a cohort that is 97% resistant genuinely
cannot yield an AUC, and that is information about the cohort.

Run `python geographic_sweep.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

MIN_ISOLATES = 30
MIN_PER_CLASS = 5


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run a geographic holdout for every eligible country.")
    parser.add_argument("--data", required=True, type=Path,
                        help="Prepared matrix containing a geographic column.")
    parser.add_argument("--config", required=True, type=Path,
                        help="Base amr.py config; split settings are overridden per country.")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Directory for the per-country runs and the summary.")
    parser.add_argument("--geographic-column", default="country")
    parser.add_argument("--min-isolates", type=int, default=MIN_ISOLATES,
                        help="Skip countries with fewer isolates than this.")
    parser.add_argument("--min-per-class", type=int, default=MIN_PER_CLASS,
                        help="Report a country as degenerate below this many of either class.")
    parser.add_argument("--python", default=sys.executable,
                        help="Interpreter used to invoke amr.py.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.out_dir.exists():
        raise SystemExit(f"{args.out_dir} already exists. Choose a new output directory.")
    if not args.data.is_file():
        raise SystemExit(f"Matrix not found: {args.data}")

    base_config = json.loads(args.config.read_text(encoding="utf-8"))
    column = args.geographic_column
    if column not in base_config.get("geographic_columns", []):
        raise SystemExit(f"'{column}' is not declared in the config's geographic_columns.")

    frame = pd.read_csv(args.data, usecols=["isolate_id", column, "resistant"],
                        dtype=str, keep_default_na=False)
    counts = (frame[frame[column] != ""]
              .groupby(column)["resistant"]
              .agg(n="count", resistant=lambda s: int((s == "1").sum())))
    counts["susceptible"] = counts["n"] - counts["resistant"]

    args.out_dir.mkdir(parents=True)
    rows: list[dict[str, object]] = []

    for country, record in counts.sort_values("n", ascending=False).iterrows():
        n = int(record["n"])
        resistant = int(record["resistant"])
        susceptible = int(record["susceptible"])
        entry: dict[str, object] = {
            "country": country, "held_out_isolates": n,
            "resistant": resistant, "susceptible": susceptible,
            "resistant_fraction": round(resistant / n, 4) if n else None,
        }

        if n < args.min_isolates:
            entry["status"] = f"skipped: fewer than {args.min_isolates} isolates"
            rows.append(entry)
            continue
        if min(resistant, susceptible) < args.min_per_class:
            # Kept in the table on purpose: a near-single-class cohort is a real
            # property of the data, not a run that merely failed.
            entry["status"] = (f"degenerate: only {min(resistant, susceptible)} of the "
                               "minority class, AUC would be uninformative")
            rows.append(entry)
            continue

        safe = "".join(ch if ch.isalnum() else "_" for ch in str(country))
        run_dir = args.out_dir / f"holdout_{safe}"
        config = dict(base_config)
        config["split_method"] = "geography_holdout"
        config["holdout_column"] = column
        config["holdout_values"] = [country]
        config["permutation_repeats"] = 0
        config_path = args.out_dir / f"config_{safe}.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

        print(f"Holding out {country} ({n} isolates, {resistant} R / {susceptible} S)...",
              file=sys.stderr, flush=True)
        completed = subprocess.run(
            [args.python, "amr.py", "train", "--data", str(args.data),
             "--config", str(config_path), "--out", str(run_dir)],
            capture_output=True, text=True)

        if completed.returncode != 0:
            message = (completed.stderr or completed.stdout or "").strip().splitlines()
            entry["status"] = "failed: " + (message[-1] if message else "unknown error")
            rows.append(entry)
            continue

        report = json.loads((run_dir / "report.json").read_text(encoding="utf-8"))
        models = report.get("models", {})
        entry["status"] = "ok"
        entry["run"] = str(run_dir)
        for variant in ("genomic_only", "full"):
            metrics = models.get(variant, {}).get("metrics")
            if metrics:
                entry[f"{variant}_auc"] = metrics.get("roc_auc")
                entry[f"{variant}_sensitivity"] = metrics.get("sensitivity")
                entry[f"{variant}_specificity"] = metrics.get("specificity")
        rows.append(entry)

    summary = pd.DataFrame(rows)
    summary_path = args.out_dir / "geographic_sweep_summary.csv"
    summary.to_csv(summary_path, index=False)
    (args.out_dir / "geographic_sweep.json").write_text(
        json.dumps({"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                    "data": str(args.data), "geographic_column": column,
                    "results": rows}, indent=2) + "\n", encoding="utf-8")

    print()
    header = f"{'country':18s} {'n':>5s} {'R':>5s} {'S':>5s} {'baseline':>9s} {'full':>8s}  status"
    print(header)
    print("-" * len(header))
    for row in rows:
        baseline = row.get("genomic_only_auc")
        full = row.get("full_auc")
        baseline_text = f"{baseline:9.3f}" if baseline is not None else f"{'-':>9s}"
        full_text = f"{full:8.3f}" if full is not None else f"{'-':>8s}"
        print(f"{str(row['country'])[:18]:18s} {row['held_out_isolates']:5d} "
              f"{row['resistant']:5d} {row['susceptible']:5d} "
              f"{baseline_text} {full_text}  {row['status']}")

    print(f"\nWrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
