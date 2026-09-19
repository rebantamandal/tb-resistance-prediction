"""Re-point a prepared matrix's grouping column without rebuilding it.

`amr.py` keeps every row sharing a `patient_id` on one side of the train/test
split. Which column fills `patient_id` is the whole experiment: group by patient
and the model may still recognise the collection a sample came from; group by
collection and it cannot; group by genetic lineage and it cannot rely on close
relatives being in training.

Building a CRyPTIC matrix means a single pass over a 1.4 GB compressed table, so
rebuilding it once per grouping wastes most of an hour for a change to one
column. The matrix already carries `site`, `subject` and `lineage`, so this
copies it and swaps `patient_id` for one of them.

    python regroup_matrix.py --source data/cryptic_rif_subject --by lineage

Run `python regroup_matrix.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd

GROUPABLE = ["subject", "site", "lineage"]
MATRIX = "isolate_variant_matrix.csv"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Copy a prepared matrix with a different grouping column.")
    parser.add_argument("--source", required=True, type=Path,
                        help="Directory holding a prepared matrix.")
    parser.add_argument("--by", required=True, choices=GROUPABLE,
                        help="Column to group by.")
    parser.add_argument("--out-dir", type=Path,
                        help="Destination; defaults to the source name with the "
                             "grouping swapped.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    matrix_path = args.source / MATRIX
    if not matrix_path.is_file():
        raise SystemExit(f"No matrix at {matrix_path}.")

    out = args.out_dir
    if out is None:
        name = args.source.name
        for candidate in GROUPABLE:
            if name.endswith("_" + candidate):
                name = name[: -len(candidate) - 1]
                break
        out = args.source.parent / f"{name}_{args.by}"
    if out.exists():
        raise SystemExit(f"{out} already exists. Choose another directory.")

    frame = pd.read_csv(matrix_path, low_memory=False,
                        dtype={c: str for c in GROUPABLE})
    if args.by not in frame.columns:
        raise SystemExit(f"The matrix has no '{args.by}' column. "
                         f"It has: {[c for c in GROUPABLE if c in frame.columns]}")
    if "patient_id" not in frame.columns:
        raise SystemExit("The matrix has no 'patient_id' column to replace.")

    source = frame[args.by].fillna("").astype(str)
    # A blank grouping value must not pool every unknown row into one enormous
    # group, which would quietly undo the split it is meant to enforce.
    blank = source.str.strip() == ""
    frame["patient_id"] = source.where(~blank, f"no{args.by}_" + frame["isolate_id"])

    out.mkdir(parents=True)
    frame.to_csv(out / MATRIX, index=False)
    for extra in ("config.generated.json", "feature_dictionary.csv"):
        if (args.source / extra).is_file():
            shutil.copy(args.source / extra, out / extra)

    report_path = args.source / "prep_report.json"
    report = json.loads(report_path.read_text(encoding="utf-8")) if report_path.is_file() else {}
    report["grouped_by"] = args.by
    report["distinct_groups"] = int(frame["patient_id"].nunique())
    report["regrouped_from"] = str(args.source)
    report["rows_with_no_grouping_value"] = int(blank.sum())
    (out / "prep_report.json").write_text(json.dumps(report, indent=2) + "\n",
                                          encoding="utf-8")

    print(f"Wrote {out / MATRIX}")
    print(f"  grouped by {args.by}: {report['distinct_groups']:,} groups "
          f"over {len(frame):,} samples")
    if blank.any():
        print(f"  {int(blank.sum()):,} samples have no {args.by} recorded and were "
              "each given their own group")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
