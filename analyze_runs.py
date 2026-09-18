"""Attach uncertainty to the metrics that `amr.py` reports as point estimates.

`amr.py` deliberately reports a single held-out number with no confidence
interval, and its own documentation says the baseline-versus-full difference is
descriptive only. That is the right default, but it makes results easy to
over-read: a difference of +0.017 AUC measured on 27 positive cases is not
distinguishable from zero, and saying so requires an interval.

This script reads the `heldout_predictions.csv` a run already wrote and
bootstraps it. Because both models scored the *same* held-out rows, the
difference is bootstrapped as a paired quantity: each resample draws a set of
isolates and recomputes both models on exactly those isolates. That is what
makes the interval on the difference meaningful rather than a subtraction of two
independent intervals.

No model is refitted here. Bootstrapping resampled *rows* measures the precision
of the estimate on this test set; it does not capture the variability that comes
from refitting on a different training split, which is larger. Treat these
intervals as a floor on the uncertainty, not a ceiling.

Run `python analyze_runs.py --help` for options.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

SCORE_SUFFIX = "_score_uncalibrated"
ACTUAL_COLUMN = "actual_resistant"


def find_models(frame: pd.DataFrame) -> list[str]:
    return [c[:-len(SCORE_SUFFIX)] for c in frame.columns if c.endswith(SCORE_SUFFIX)]


def bootstrap_indices(n: int, resamples: int, rng: np.random.Generator) -> np.ndarray:
    return rng.integers(0, n, size=(resamples, n))


def summarise(values: np.ndarray, level: float) -> dict[str, float | None]:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"point": None, "low": None, "high": None}
    tail = (1.0 - level) / 2.0
    return {
        "median": float(np.median(values)),
        "low": float(np.quantile(values, tail)),
        "high": float(np.quantile(values, 1.0 - tail)),
    }


def analyse_run(run_dir: Path, resamples: int, level: float,
                seed: int) -> dict[str, object]:
    predictions_path = run_dir / "heldout_predictions.csv"
    if not predictions_path.is_file():
        raise SystemExit(f"{predictions_path} not found; is this an amr.py run directory?")

    frame = pd.read_csv(predictions_path)
    if ACTUAL_COLUMN not in frame.columns:
        raise SystemExit(f"{predictions_path} has no '{ACTUAL_COLUMN}' column.")

    actual = frame[ACTUAL_COLUMN].to_numpy(dtype=float)
    models = find_models(frame)
    if not models:
        raise SystemExit(f"{predictions_path} contains no score columns.")

    scores = {m: frame[f"{m}{SCORE_SUFFIX}"].to_numpy(dtype=float) for m in models}
    n = len(actual)
    positives = int(actual.sum())

    rng = np.random.default_rng(seed)
    draws = bootstrap_indices(n, resamples, rng)

    point = {m: {"roc_auc": float(roc_auc_score(actual, scores[m])),
                 "average_precision": float(average_precision_score(actual, scores[m]))}
             for m in models}

    boot_auc: dict[str, list[float]] = {m: [] for m in models}
    differences: dict[str, list[float]] = {}
    baseline = "genomic_only" if "genomic_only" in models else models[0]
    others = [m for m in models if m != baseline]
    for other in others:
        differences[other] = []

    for draw in draws:
        resampled_actual = actual[draw]
        # A resample containing only one class has no defined AUC; skip it
        # rather than substituting a value, and report how many were skipped.
        if resampled_actual.min() == resampled_actual.max():
            continue
        per_model = {}
        for model in models:
            value = roc_auc_score(resampled_actual, scores[model][draw])
            per_model[model] = value
            boot_auc[model].append(value)
        for other in others:
            differences[other].append(per_model[other] - per_model[baseline])

    usable = len(boot_auc[baseline])
    result: dict[str, object] = {
        "run": str(run_dir),
        "held_out_isolates": n,
        "held_out_resistant": positives,
        "held_out_susceptible": n - positives,
        "bootstrap_resamples_requested": resamples,
        "bootstrap_resamples_usable": usable,
        "confidence_level": level,
        "models": {},
        "comparisons": {},
    }
    for model in models:
        result["models"][model] = {
            "roc_auc_point": point[model]["roc_auc"],
            "average_precision_point": point[model]["average_precision"],
            "roc_auc_ci": summarise(np.array(boot_auc[model]), level),
        }
    for other in others:
        values = np.array(differences[other])
        summary = summarise(values, level)
        summary["point"] = point[other]["roc_auc"] - point[baseline]["roc_auc"]
        # The share of resamples favouring the richer model. This is a direct
        # statement of how consistent the improvement is, not a p-value.
        summary["fraction_of_resamples_favouring_" + other] = (
            float((values > 0).mean()) if values.size else None)
        summary["interval_excludes_zero"] = bool(
            summary["low"] is not None and (summary["low"] > 0 or summary["high"] < 0))
        result["comparisons"][f"{other} minus {baseline}"] = summary
    return result


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bootstrap confidence intervals for saved amr.py runs.")
    parser.add_argument("--runs", required=True, nargs="+", type=Path,
                        help="One or more amr.py run directories.")
    parser.add_argument("--out", type=Path, help="Optional JSON summary path.")
    parser.add_argument("--resamples", type=int, default=2000)
    parser.add_argument("--confidence", type=float, default=0.95)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)
    if not 100 <= args.resamples <= 100_000:
        parser.error("--resamples must be between 100 and 100000.")
    if not 0.5 < args.confidence < 1.0:
        parser.error("--confidence must be strictly between 0.5 and 1.")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    results = [analyse_run(run, args.resamples, args.confidence, args.seed)
               for run in args.runs]

    percent = int(args.confidence * 100)
    print(f"{'run':38s} {'model':14s} {'n':>5s} {'R':>4s} "
          f"{'AUC':>7s}  {percent}% interval")
    print("-" * 88)
    for result in results:
        for model, stats in result["models"].items():
            interval = stats["roc_auc_ci"]
            span = (f"[{interval['low']:.3f}, {interval['high']:.3f}]"
                    if interval["low"] is not None else "undefined")
            print(f"{Path(result['run']).name:38s} {model:14s} "
                  f"{result['held_out_isolates']:5d} {result['held_out_resistant']:4d} "
                  f"{stats['roc_auc_point']:7.3f}  {span}")

    comparisons = [(r, k, v) for r in results for k, v in r["comparisons"].items()]
    if comparisons:
        print()
        print(f"{'run':38s} {'comparison':26s} {'diff':>7s}  {percent}% interval")
        print("-" * 88)
        for result, name, stats in comparisons:
            span = (f"[{stats['low']:+.3f}, {stats['high']:+.3f}]"
                    if stats["low"] is not None else "undefined")
            verdict = "excludes 0" if stats["interval_excludes_zero"] else "includes 0"
            print(f"{Path(result['run']).name:38s} {name:26s} "
                  f"{stats['point']:+7.3f}  {span}  {verdict}")

    if args.out:
        if args.out.exists():
            raise SystemExit(f"{args.out} already exists. Choose a new filename.")
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
        print(f"\nWrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
