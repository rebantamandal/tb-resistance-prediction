"""Run the complete evaluation ladder for one or more drugs.

For each drug this performs, in order:

  1. an isolate-grouped genomic-only baseline      (the permissive protocol)
  2. a study-grouped baseline and full model       (the strict protocol)
  3. a geographic holdout for every eligible country

The point of running all three is that only their comparison is informative. A
single number from step 1 is the figure usually published; steps 2 and 3 say how
much of it survives a harder question.

`host` is deliberately excluded from the clinical block throughout. Its
missingness identifies one high-resistance cohort, so including it lets the model
recognise the study rather than the biology (see RESULTS.md §3).

Run `python run_full_study.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PYTHON = sys.executable
CLINICAL = ["infection_site"]


def run(command: list[str], label: str) -> bool:
    print(f"\n=== {label} ===", flush=True)
    completed = subprocess.run(command, capture_output=True, text=True)
    output = (completed.stdout or "") + (completed.stderr or "")
    for line in output.strip().splitlines()[-6:]:
        print("   ", line, flush=True)
    if completed.returncode != 0:
        print(f"    FAILED ({completed.returncode})", flush=True)
    return completed.returncode == 0


def write_run_config(prepared: Path, enable_full: bool) -> Path:
    """Take the generated config and drop `host` from the clinical block."""
    config = json.loads((prepared / "config.generated.json").read_text(encoding="utf-8"))
    if enable_full:
        config["clinical_categorical_columns"] = [
            c for c in CLINICAL if c in config.get("clinical_categorical_columns", [])]
        config["compare_full"] = bool(config["clinical_categorical_columns"]
                                      and config.get("geographic_columns"))
    path = prepared / "config.run.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the full evaluation ladder.")
    parser.add_argument("--drugs", required=True, nargs="+",
                        help="Drug column names in the phenotype table, e.g. ETHAMBUTOL.")
    parser.add_argument("--variants", type=Path, default=Path("all_variants.csv"))
    parser.add_argument("--phenotypes", type=Path, default=Path("data/phenotypes.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/sample_metadata.csv"))
    parser.add_argument("--min-prevalence", type=int, default=20)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--models-dir", type=Path, default=Path("models"))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    for drug in args.drugs:
        slug = drug.lower()[:3]
        pretty = drug.capitalize()
        base = ["--variants", str(args.variants), "--min-prevalence", str(args.min_prevalence),
                "--phenotypes", str(args.phenotypes), "--drug-column", drug,
                "--antibiotic", pretty, "--quiet"]

        # 1. permissive: grouped by isolate, genomic only
        prepared = args.data_dir / f"{slug}_baseline"
        if not prepared.exists():
            run([PYTHON, "prepare_variants.py", *base, "--out-dir", str(prepared)],
                f"{pretty}: prepare isolate-grouped matrix")
        out = args.models_dir / f"{slug}_genomic_baseline"
        if prepared.exists() and not out.exists():
            run([PYTHON, "amr.py", "train", "--data",
                 str(prepared / "isolate_variant_matrix.csv"),
                 "--config", str(prepared / "config.generated.json"), "--out", str(out)],
                f"{pretty}: train isolate-grouped baseline")

        # 2. strict: grouped by study, baseline versus full
        prepared = args.data_dir / f"{slug}_bystudy"
        if not prepared.exists():
            run([PYTHON, "prepare_variants.py", *base, "--out-dir", str(prepared),
                 "--metadata", str(args.metadata), "--group-by", "study"],
                f"{pretty}: prepare study-grouped matrix")
        out = args.models_dir / f"{slug}_bystudy"
        if prepared.exists() and not out.exists():
            config = write_run_config(prepared, enable_full=True)
            run([PYTHON, "amr.py", "train", "--data",
                 str(prepared / "isolate_variant_matrix.csv"),
                 "--config", str(config), "--out", str(out)],
                f"{pretty}: train study-grouped baseline and full")

        # 3. geographic holdouts, one country at a time
        prepared = args.data_dir / f"{slug}_geo"
        if not prepared.exists():
            run([PYTHON, "prepare_variants.py", *base, "--out-dir", str(prepared),
                 "--metadata", str(args.metadata), "--group-by", "study",
                 "--require-geography"],
                f"{pretty}: prepare country-complete matrix")
        out = args.models_dir / f"geo_sweep_{slug}"
        if prepared.exists() and not out.exists():
            config = write_run_config(prepared, enable_full=True)
            run([PYTHON, "geographic_sweep.py", "--data",
                 str(prepared / "isolate_variant_matrix.csv"),
                 "--config", str(config), "--out-dir", str(out)],
                f"{pretty}: geographic holdout sweep")

    print("\nAll requested drugs complete.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
