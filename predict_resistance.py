"""Predict antibiotic resistance for new isolates from raw variant calls.

This is the end-to-end tool: give it a variant file in the same shape as the
source data and it returns, for every isolate, a resistant/susceptible call per
drug with a score.

    python predict_resistance.py --variants new_isolates.csv --out predictions.csv

The input is the same five columns the training data uses:

    SAMPLE,CHROM,POS,REF,ALT
    ERR038266,NC_000962.3,1849,C,A

One row per variant call per isolate, any number of isolates. The script builds
exactly the feature columns each saved model was trained on — same positions,
same alleles, same order — then runs every model registered in the model
directory.

A variant the model has never seen is ignored, because the model has no column
for it. An isolate carrying no known feature at all is flagged rather than
silently scored, since an all-zero row is indistinguishable from a fully
susceptible one and the two mean very different things.

RESEARCH USE ONLY. Scores are uncalibrated and the models are not clinically
validated. A score of 0.80 is not an 80% probability that this patient's
infection is resistant. Do not use this to start, stop or change treatment.

Run `python predict_resistance.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

CHUNK_ROWS = 2_000_000
VARIANT_COLUMNS = ["SAMPLE", "CHROM", "POS", "REF", "ALT"]
NOTICE = ("RESEARCH ONLY. Uncalibrated scores, not clinically validated. "
          "Not a patient's probability of resistance. Do not use for treatment decisions.")

# Each entry: the drug, the run directory holding its model, and the held-out
# performance actually measured for it. The optimistic and realistic figures are
# both carried so the output can state what the number is worth.
DEFAULT_REGISTRY = {
    "Rifampicin": {
        "run": "models/rif_genomic_baseline",
        "auc_random_split": 0.958,
        "auc_study_grouped": 0.908,
        "auc_unseen_country_range": "0.69-0.76",
    },
    "Isoniazid": {
        "run": "models/iso_genomic_baseline",
        "auc_random_split": 0.947,
        "auc_study_grouped": 0.855,
        "auc_unseen_country_range": "0.77-0.94",
    },
    "Ethambutol": {
        "run": "models/eth_genomic_baseline",
        "auc_random_split": 0.913,
        "auc_study_grouped": 0.765,
        "auc_unseen_country_range": "0.63-0.84",
    },
    "Pyrazinamide": {
        "run": "models/pyr_genomic_baseline",
        "auc_random_split": 0.911,
        "auc_study_grouped": 0.782,
        "auc_unseen_country_range": "0.63-0.84",
    },
}


def variant_key(frame: pd.DataFrame) -> pd.Series:
    return (frame["POS"].astype(str) + "_" + frame["REF"].astype(str)
            + "_" + frame["ALT"].astype(str))


def load_bundle(model_path: Path) -> dict:
    """Load a saved amr.py model bundle: the fitted pipeline plus its metadata."""
    bundle = joblib.load(model_path)
    if not isinstance(bundle, dict) or "pipeline" not in bundle:
        raise SystemExit(f"{model_path} is not an amr.py model bundle.")
    return bundle


def bundle_feature_columns(bundle: dict, model_path: Path) -> list[str]:
    """The exact genomic columns this model expects, in training order.

    Taken from the bundle's own input schema rather than a config file beside
    it, so the columns always come from the object actually being scored.
    """
    schema = bundle.get("input_schema") or []
    columns = [entry["name"] for entry in schema if entry.get("kind") == "genomic"]
    if not columns:
        columns = list(bundle.get("config", {}).get("genomic_columns") or [])
    if not columns:
        raise SystemExit(f"{model_path} declares no genomic feature columns.")
    return columns


def check_versions(bundle: dict, model_path: Path, drug: str) -> list[str]:
    """Warn when the current environment differs from the one that fitted the model.

    A pickled estimator loaded under a different scikit-learn can score
    differently or fail subtly, so a mismatch is surfaced rather than ignored.
    """
    import sklearn
    trained = bundle.get("versions") or {}
    current = {"scikit-learn": sklearn.__version__, "pandas": pd.__version__,
               "numpy": np.__version__}
    return [f"{drug}: trained with {name} {trained[name]}, running {value}"
            for name, value in current.items()
            if name in trained and trained[name] != value]


def column_to_key(columns: list[str], dictionary: pd.DataFrame) -> dict[str, str]:
    """Map each model column back to the variant key that switches it on.

    The feature dictionary written alongside the prepared matrix is the
    authority here: column names encode gene and position but shorten long
    indel alleles, so they cannot be parsed back into a key reliably.
    """
    lookup = {}
    for feature, position, ref, alt in zip(dictionary["feature"], dictionary["position"],
                                           dictionary["ref"], dictionary["alt"]):
        lookup[str(feature)] = f"{position}_{ref}_{alt}"
    missing = [c for c in columns if c not in lookup]
    if missing:
        raise SystemExit(
            f"{len(missing)} model features are absent from the feature dictionary "
            f"(first: {missing[0]}). Point --feature-dictionary at the dictionary "
            "written by the same prepare_variants.py run that produced these models.")
    return {c: lookup[c] for c in columns}


def build_matrix(variants_path: Path, key_of_column: dict[str, str],
                 progress: bool) -> tuple[pd.DataFrame, pd.Series]:
    """Stream the variant file and switch on the columns the models know about."""
    columns = list(key_of_column)
    index_of_key: dict[str, list[int]] = {}
    for position, column in enumerate(columns):
        index_of_key.setdefault(key_of_column[column], []).append(position)

    samples: list[str] = []
    sample_index: dict[str, int] = {}
    rows: list[np.ndarray] = []
    calls_seen: dict[str, int] = {}
    rows_read = 0

    reader = pd.read_csv(variants_path, chunksize=CHUNK_ROWS, dtype=str,
                         keep_default_na=False, na_filter=False)
    for chunk in reader:
        missing = [c for c in VARIANT_COLUMNS if c not in chunk.columns]
        if missing:
            raise SystemExit(f"{variants_path} is missing columns: {missing}. "
                             f"Expected {VARIANT_COLUMNS}.")
        rows_read += len(chunk)
        keys = variant_key(chunk)
        for sample, key in zip(chunk["SAMPLE"].to_numpy(), keys.to_numpy()):
            position = sample_index.get(sample)
            if position is None:
                position = len(samples)
                sample_index[sample] = position
                samples.append(sample)
                rows.append(np.zeros(len(columns), dtype=np.uint8))
                calls_seen[sample] = 0
            calls_seen[sample] += 1
            for column_position in index_of_key.get(key, ()):
                rows[position][column_position] = 1
        if progress:
            print(f"  read {rows_read:,} variant rows, {len(samples):,} isolates",
                  file=sys.stderr, flush=True)

    if not samples:
        raise SystemExit(f"No variant rows found in {variants_path}.")
    matrix = pd.DataFrame(np.vstack(rows), columns=columns)
    matrix.insert(0, "isolate_id", samples)
    return matrix, pd.Series([calls_seen[s] for s in samples], index=range(len(samples)))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Predict antibiotic resistance from raw variant calls.")
    parser.add_argument("--variants", required=True, type=Path,
                        help="CSV with SAMPLE,CHROM,POS,REF,ALT for the isolates to score.")
    parser.add_argument("--out", required=True, type=Path,
                        help="Output CSV; refuses to overwrite an existing file.")
    parser.add_argument("--feature-dictionary", type=Path,
                        default=Path("data/variant_panel/feature_dictionary.csv"),
                        help="Dictionary mapping feature columns to variant positions.")
    parser.add_argument("--registry", type=Path,
                        help="Optional JSON overriding which drugs and runs to use.")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="Score at or above which an isolate is called resistant. "
                             "Fixed at 0.50 during evaluation; changing it changes "
                             "the sensitivity/specificity trade-off, not the model.")
    parser.add_argument("--trust-local-models", action="store_true",
                        help="Required. Model files execute code when loaded, so confirm "
                             "these were produced locally by you.")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    if not 0.0 < args.threshold < 1.0:
        parser.error("--threshold must be strictly between 0 and 1.")
    if not args.trust_local_models:
        parser.error("Pass --trust-local-models to confirm the model files are yours. "
                     "Never load a model file from an untrusted source.")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    progress = not args.quiet
    print(NOTICE)

    if args.out.exists():
        raise SystemExit(f"{args.out} already exists. Choose a new filename.")
    if not args.variants.is_file():
        raise SystemExit(f"Variant file not found: {args.variants}")
    if not args.feature_dictionary.is_file():
        raise SystemExit(f"Feature dictionary not found: {args.feature_dictionary}")

    registry = (json.loads(args.registry.read_text(encoding="utf-8"))
                if args.registry else DEFAULT_REGISTRY)

    available = {}
    for drug, entry in registry.items():
        run_dir = Path(entry["run"])
        model_path = run_dir / "genomic_only.joblib"
        if model_path.is_file():
            available[drug] = {**entry, "run_dir": run_dir, "model_path": model_path}
        elif progress:
            print(f"  skipping {drug}: {model_path} not found", file=sys.stderr)
    if not available:
        raise SystemExit("No trained models found. Train at least one drug first.")

    # Load every bundle up front so a broken or mismatched model fails before
    # the expensive pass over the variant file rather than after it.
    version_warnings: list[str] = []
    for drug, entry in available.items():
        entry["bundle"] = load_bundle(entry["model_path"])
        entry["columns"] = bundle_feature_columns(entry["bundle"], entry["model_path"])
        version_warnings.extend(check_versions(entry["bundle"], entry["model_path"], drug))

    # Every registered model must expect the same feature set, or one matrix
    # cannot serve them all.
    reference_drug = next(iter(available))
    reference_columns = available[reference_drug]["columns"]
    for drug, entry in available.items():
        if entry["columns"] != reference_columns:
            raise SystemExit(
                f"{drug} expects different features from {reference_drug}. "
                "Models scored together must come from the same feature panel. "
                "Re-run prepare_variants.py and retrain so both share one panel.")

    if version_warnings:
        print("\nEnvironment differs from the one these models were fitted in:")
        for warning in version_warnings:
            print(f"  {warning}")
        print("  Scores may differ from the evaluated models. Retrain, or pin the "
              "versions in requirements.txt.")

    dictionary = pd.read_csv(args.feature_dictionary, dtype=str)
    key_of_column = column_to_key(reference_columns, dictionary)

    if progress:
        print(f"Scoring with {len(available)} model(s): {', '.join(sorted(available))}",
              file=sys.stderr)
    matrix, calls_per_isolate = build_matrix(args.variants, key_of_column, progress)

    features = matrix[reference_columns]
    known_features_present = features.sum(axis=1)

    result = pd.DataFrame({
        "isolate_id": matrix["isolate_id"],
        "variant_calls_supplied": calls_per_isolate.to_numpy(),
        "known_features_present": known_features_present.to_numpy(),
    })

    for drug in sorted(available):
        entry = available[drug]
        scores = entry["bundle"]["pipeline"].predict_proba(features)[:, 1]
        result[f"{drug}_score"] = np.round(scores, 4)
        result[f"{drug}_call"] = np.where(scores >= args.threshold,
                                          "Resistant", "Susceptible")

    # An isolate with no recognised feature scores like a susceptible one for
    # purely structural reasons. Say so rather than letting the call stand alone.
    result["warning"] = np.where(
        result["known_features_present"] == 0,
        "No known resistance feature detected; an all-zero profile is not evidence "
        "of susceptibility. Check that coordinates match the H37Rv reference.", "")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "notice": NOTICE,
        "variant_file": str(args.variants),
        "isolates_scored": len(result),
        "threshold": args.threshold,
        "feature_count": len(reference_columns),
        "isolates_with_no_known_feature": int((result["known_features_present"] == 0).sum()),
        "models": {drug: {
            "run": str(entry["run_dir"]),
            "held_out_auc_random_split": entry.get("auc_random_split"),
            "held_out_auc_study_grouped": entry.get("auc_study_grouped"),
            "auc_on_an_unseen_country": entry.get("auc_unseen_country_range"),
        } for drug, entry in sorted(available.items())},
        "how_to_read_the_performance_figures": (
            "The random-split figure is the one usually quoted and is the most "
            "optimistic. The study-grouped figure is what to expect on isolates "
            "from a cohort the model did not train on, and is the figure to plan "
            "around. The unseen-country entry is a RANGE across several countries, "
            "not a floor: it is noisy and not uniformly worse than the "
            "study-grouped figure. Rifampicin drops on every country tested, while "
            "isoniazid scores above its study-grouped value on two of four. Treat "
            "any single country figure as weak evidence; see RESULTS.md section 2."),
    }
    report_path = args.out.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"\nWrote {args.out} ({len(result):,} isolates scored).")
    for drug in sorted(available):
        calls = result[f"{drug}_call"].value_counts()
        print(f"  {drug:14s} {calls.get('Resistant', 0):5d} resistant / "
              f"{calls.get('Susceptible', 0):5d} susceptible")
    blank = int((result["known_features_present"] == 0).sum())
    if blank:
        print(f"\n  {blank} isolate(s) carried no known resistance feature and are flagged.")
    print(f"Report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
