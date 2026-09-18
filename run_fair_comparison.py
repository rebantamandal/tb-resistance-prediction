"""The fair test of the project's hypothesis: do clinical and geographic features help?

Two earlier comparisons each had a defect that made them uninformative:

  1. Complete metadata but grouped by isolate. Isolates from one study sat on
     both sides of the split, so country could act as a lookup for that study's
     resistance rate. The full model looked better than it was.

  2. Grouped by study but with metadata missing for most rows. Country was blank
     for 78% of isolates, and blankness itself predicted resistance (61.8%
     resistant when missing, 30.6% when present). The full model was handed a
     feature whose absence encoded cohort, and looked worse than it was.

This script runs the cell neither of those covered: **complete metadata and
grouped by study**. Every isolate has a country and an infection site, so there
is no missingness signal, and no study spans the split, so there is no cohort
lookup. Whatever difference remains is the honest answer.

Permutation importance is enabled so the added features can be read individually
rather than only as a block.

Run `python run_fair_comparison.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PYTHON = sys.executable
CLINICAL = ["infection_site"]
GEOGRAPHIC = ["country", "city_or_region"]


def run(command: list[str], label: str) -> bool:
    print(f"\n=== {label} ===", flush=True)
    completed = subprocess.run(command, capture_output=True, text=True)
    text = (completed.stdout or "") + (completed.stderr or "")
    for line in text.strip().splitlines()[-4:]:
        print("   ", line, flush=True)
    if completed.returncode != 0:
        print(f"    FAILED ({completed.returncode})", flush=True)
    return completed.returncode == 0


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Baseline vs full on complete metadata, grouped by study.")
    parser.add_argument("--drugs", required=True, nargs="+")
    parser.add_argument("--variants", type=Path, default=Path("all_variants.csv"))
    parser.add_argument("--phenotypes", type=Path, default=Path("data/phenotypes.csv"))
    parser.add_argument("--metadata", type=Path, default=Path("data/sample_metadata.csv"))
    parser.add_argument("--min-prevalence", type=int, default=20)
    parser.add_argument("--permutation-repeats", type=int, default=5)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    for drug in args.drugs:
        slug = drug.lower()[:3]
        prepared = Path("data") / f"{slug}_fair"
        out = Path("models") / f"{slug}_fair"

        if not prepared.exists():
            ok = run([PYTHON, "prepare_variants.py",
                      "--variants", str(args.variants),
                      "--min-prevalence", str(args.min_prevalence),
                      "--phenotypes", str(args.phenotypes),
                      "--drug-column", drug,
                      "--antibiotic", drug.capitalize(),
                      "--metadata", str(args.metadata),
                      "--group-by", "study",
                      "--require-metadata",
                      "--out-dir", str(prepared), "--quiet"],
                     f"{drug}: prepare complete-metadata, study-grouped matrix")
            if not ok:
                continue

        config = json.loads((prepared / "config.generated.json").read_text(encoding="utf-8"))
        # host stays out: its missingness identifies one high-resistance cohort.
        config["clinical_categorical_columns"] = [
            c for c in CLINICAL if c in config.get("clinical_categorical_columns", [])]
        config["geographic_columns"] = [
            c for c in GEOGRAPHIC if c in config.get("geographic_columns", [])]
        config["compare_full"] = bool(config["clinical_categorical_columns"]
                                      and config["geographic_columns"])
        config["permutation_repeats"] = args.permutation_repeats
        config_path = prepared / "config.fair.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

        if not config["compare_full"]:
            print(f"    {drug}: no usable clinical/geographic columns, skipping")
            continue

        if not out.exists():
            run([PYTHON, "amr.py", "train",
                 "--data", str(prepared / "isolate_variant_matrix.csv"),
                 "--config", str(config_path), "--out", str(out)],
                f"{drug}: train baseline vs full on identical rows")

    print("\nDone.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
