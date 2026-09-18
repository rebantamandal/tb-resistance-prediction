"""Research-only antibiotic resistance classifier. No pretrained model is included.

This module consumes a prepared tabular feature matrix, not raw sequencing files.
It fits one organism-antibiotic pair at a time and compares feature sets on an
identical group-disjoint holdout. Run `python amr.py --help` for the CLI.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.metrics import (accuracy_score, average_precision_score,
                             confusion_matrix, f1_score, roc_auc_score)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

NOTICE = ("RESEARCH ONLY. Not clinically validated. Scores are uncalibrated, "
          "not a patient's proven probability of resistance. Do not use for treatment decisions.")
FORMAT_VERSION = 1


@dataclass
class Config:
    organism: str
    antibiotic: str
    genomic_columns: list[str]
    organism_column: str = "organism"
    antibiotic_column: str = "antibiotic"
    target_column: str = "resistant"
    sample_column: str = "isolate_id"
    group_column: str = "patient_id"
    clinical_numeric_columns: list[str] = field(default_factory=list)
    clinical_categorical_columns: list[str] = field(default_factory=list)
    geographic_columns: list[str] = field(default_factory=list)
    compare_full: bool = True
    split_method: str = "group_random"
    holdout_column: str = ""
    holdout_values: list[str] = field(default_factory=list)
    test_size: float = 0.2
    threshold: float = 0.5
    n_estimators: int = 300
    random_state: int = 42
    permutation_repeats: int = 0

    def validate(self) -> None:
        for value in (self.organism, self.antibiotic, self.organism_column,
                      self.antibiotic_column, self.target_column,
                      self.sample_column, self.group_column):
            if not isinstance(value, str) or not value.strip():
                raise ValueError("Organism, antibiotic, target and identifier settings must be non-empty strings.")
        if self.organism.startswith("REPLACE_") or self.antibiotic.startswith("REPLACE_"):
            raise ValueError("Replace the organism and antibiotic placeholders with exact values from your data.")
        if not self.genomic_columns:
            raise ValueError("Select at least one measured genomic feature.")
        blocks = [self.genomic_columns, self.clinical_numeric_columns,
                  self.clinical_categorical_columns, self.geographic_columns]
        if any(not isinstance(b, list) for b in blocks):
            raise ValueError("Feature selections must be lists of column names.")
        features = sum(blocks, [])
        if any(not isinstance(c, str) or not c.strip() for c in features):
            raise ValueError("Feature names must be non-empty strings.")
        if len(features) != len(set(features)):
            raise ValueError("A feature cannot appear in more than one feature block.")
        protected = {self.organism_column, self.antibiotic_column, self.target_column,
                     self.sample_column, self.group_column}
        if protected.intersection(features):
            raise ValueError("Identifiers, organism, antibiotic and the target cannot be predictor columns.")
        if len({self.organism_column, self.antibiotic_column, self.target_column}) != 3:
            raise ValueError("Organism, antibiotic and target must be different columns.")
        if {self.sample_column, self.group_column}.intersection(
                {self.organism_column, self.antibiotic_column, self.target_column}):
            raise ValueError("Identifier columns cannot also be scope or target columns.")
        if self.compare_full and (not self.geographic_columns or not
                (self.clinical_numeric_columns or self.clinical_categorical_columns)):
            raise ValueError("Full comparison needs genuine clinical AND geographic features. Otherwise use genomic-only mode.")
        if not 0.1 <= float(self.test_size) <= 0.5:
            raise ValueError("test_size must be between 0.1 and 0.5.")
        if not 0.0 < float(self.threshold) < 1.0:
            raise ValueError("threshold must be strictly between 0 and 1.")
        if not isinstance(self.n_estimators, int) or not 10 <= self.n_estimators <= 2000:
            raise ValueError("n_estimators must be an integer between 10 and 2000.")
        if not isinstance(self.random_state, int) or not 0 <= self.random_state < 2**32:
            raise ValueError("random_state must be an integer from 0 to 2**32-1.")
        if not isinstance(self.permutation_repeats, int) or not 0 <= self.permutation_repeats <= 10:
            raise ValueError("permutation_repeats must be an integer from 0 to 10.")
        if self.split_method not in {"group_random", "geography_holdout"}:
            raise ValueError("Unknown split_method.")
        if self.split_method == "geography_holdout":
            if self.holdout_column not in self.geographic_columns or not self.holdout_values:
                raise ValueError("Choose a declared geographic column and at least one value to hold out.")

    def feature_blocks(self, variant: str) -> tuple[list[str], list[str], list[str]]:
        if variant == "genomic_only":
            return self.genomic_columns, [], []
        if variant == "full" and self.compare_full:
            return (self.genomic_columns, self.clinical_numeric_columns,
                    self.clinical_categorical_columns + self.geographic_columns)
        raise ValueError(f"Unknown or disabled model variant: {variant}")


def versions() -> dict[str, str]:
    return {"python": platform.python_version(), "scikit-learn": sklearn.__version__,
            "pandas": pd.__version__, "numpy": np.__version__, "joblib": joblib.__version__}


def read_csv(source: Any) -> pd.DataFrame:
    """Preserve string identifiers and reject ambiguous duplicate column headers."""
    if isinstance(source, (str, Path)):
        text = Path(source).read_text(encoding="utf-8-sig")
    else:
        text = source.read()
        if isinstance(text, bytes):
            text = text.decode("utf-8-sig")
    text = text.lstrip("\ufeff")
    try:
        header = next(csv.reader(io.StringIO(text)))
    except StopIteration as exc:
        raise ValueError("The CSV is empty.") from exc
    header = [x.strip() for x in header]
    if not header or any(not c for c in header) or len(header) != len(set(header)):
        raise ValueError("CSV column names must be unique and non-empty.")
    frame = pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)
    frame.columns = header
    for col in frame:
        frame[col] = frame[col].map(lambda x: x.strip() if isinstance(x, str) else x)
        frame[col] = frame[col].replace("", np.nan)
    return frame


def _require(frame: pd.DataFrame, columns: list[str]) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))
    if frame.empty:
        raise ValueError("There are no data rows. Templates contain headers only; add real observations.")


def _scoped(frame: pd.DataFrame, cfg: Config, *, training: bool) -> pd.DataFrame:
    _require(frame, [cfg.organism_column, cfg.antibiotic_column])
    selected = (frame[cfg.organism_column].eq(cfg.organism) &
                frame[cfg.antibiotic_column].eq(cfg.antibiotic))
    if not training and not selected.all():
        raise ValueError("This model only accepts its trained organism-antibiotic pair. Batch contains a different or missing pair.")
    result = frame.loc[selected].copy().reset_index(drop=True)
    if result.empty:
        raise ValueError("No rows match the selected organism and antibiotic. Check spelling and case.")
    return result


def encode_target(series: pd.Series) -> pd.Series:
    mapping = {"r": 1, "resistant": 1, "1": 1, "1.0": 1,
               "s": 0, "susceptible": 0, "0": 0, "0.0": 0}
    normalized = series.map(lambda x: str(x).strip().lower() if pd.notna(x) else "<missing>")
    invalid = sorted(set(normalized) - set(mapping))
    if invalid:
        raise ValueError("Target must be R/S, resistant/susceptible, or 1/0. Invalid values: "
                         + ", ".join(invalid[:8]) + ". Intermediate, unknown, missing and -1 labels "
                         "are NOT silently reclassified. Resolve/exclude them with a documented policy first.")
    return normalized.map(mapping).astype(int)


def prepare_features(frame: pd.DataFrame, cfg: Config, variant: str) -> pd.DataFrame:
    genes, numeric, categorical = cfg.feature_blocks(variant)
    columns = genes + numeric + categorical
    _require(frame, columns)
    result = frame[columns].copy()
    for col in genes + numeric:
        try:
            result[col] = pd.to_numeric(result[col], errors="raise").astype(float)
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{col}: use numeric values; leave genuinely unknown entries blank.") from exc
        if np.isinf(result[col]).any():
            raise ValueError(f"{col}: infinity is not a valid value.")
        if col in {"age", "hospitalization_days_before_sample"} and (result[col].dropna() < 0).any():
            raise ValueError(f"{col}: negative values are invalid; unknown values should be blank.")
        if col == "prior_antibiotic_exposure" and not result[col].dropna().isin([0.0, 1.0]).all():
            raise ValueError("prior_antibiotic_exposure: this template requires 0, 1, or blank.")
        if col in genes and not result[col].dropna().isin([0.0, 1.0]).all():
            raise ValueError(f"{col}: genomic indicators must be 0, 1, or blank. Map source-specific missing codes explicitly.")
    for col in categorical:
        result[col] = result[col].map(lambda x: str(x).strip() if pd.notna(x) and str(x).strip() else np.nan)
    if result.isna().all(axis=1).any():
        raise ValueError(f"At least one row has no observed predictors for {variant}.")
    return result


def make_pipeline(cfg: Config, variant: str) -> Pipeline:
    genes, numeric, categorical = cfg.feature_blocks(variant)
    transforms = [("genomic", SimpleImputer(strategy="most_frequent", add_indicator=True), genes)]
    if numeric:
        transforms.append(("clinical_numeric", SimpleImputer(strategy="median", add_indicator=True), numeric))
    if categorical:
        cat_pipe = Pipeline([
            ("missing", SimpleImputer(strategy="constant", fill_value="__MISSING__")),
            ("encode", OneHotEncoder(handle_unknown="ignore"))])
        transforms.append(("categorical", cat_pipe, categorical))
    return Pipeline([
        ("prepare", ColumnTransformer(transforms, remainder="drop")),
        ("forest", RandomForestClassifier(
            n_estimators=cfg.n_estimators, min_samples_leaf=2,
            class_weight="balanced_subsample", random_state=cfg.random_state,
            n_jobs=2))])


def _ratio(a: int, b: int) -> float | None:
    return float(a / b) if b else None


def metrics(y: Any, scores: Any, threshold: float) -> dict[str, Any]:
    y = np.asarray(y, dtype=int)
    scores = np.asarray(scores, dtype=float)
    predicted = (scores >= threshold).astype(int)
    tn, fp, fn, tp = [int(x) for x in confusion_matrix(y, predicted, labels=[0, 1]).ravel()]
    both = len(np.unique(y)) == 2
    return {"n": len(y), "resistant_n": int(y.sum()), "susceptible_n": int((y == 0).sum()),
            "resistance_prevalence": float(y.mean()), "accuracy": float(accuracy_score(y, predicted)),
            "sensitivity": _ratio(tp, tp + fn), "specificity": _ratio(tn, tn + fp),
            "precision": _ratio(tp, tp + fp), "f1": float(f1_score(y, predicted, zero_division=0)),
            "roc_auc": float(roc_auc_score(y, scores)) if both else None,
            "average_precision": float(average_precision_score(y, scores)) if both else None,
            "true_negative": tn, "false_positive": fp, "false_negative": fn, "true_positive": tp}


def _positive_scores(pipeline: Pipeline, X: pd.DataFrame) -> np.ndarray:
    index = list(pipeline.named_steps["forest"].classes_).index(1)
    return pipeline.predict_proba(X)[:, index]


def _input_schema(X_train: pd.DataFrame, cfg: Config, variant: str) -> list[dict[str, Any]]:
    genes, numeric, categorical = cfg.feature_blocks(variant)
    schema = []
    for col in X_train:
        item: dict[str, Any] = {"name": col, "kind": "genomic" if col in genes else
                              "numeric" if col in numeric else "categorical"}
        values = X_train[col].dropna()
        if col in categorical:
            item["values"] = sorted(values.unique().tolist())
        else:
            item["min"] = float(values.min())
            item["max"] = float(values.max())
        schema.append(item)
    return schema


def split_rows(frame: pd.DataFrame, cfg: Config, y: pd.Series) -> tuple[np.ndarray, np.ndarray, int]:
    groups = frame[cfg.group_column].astype(str)
    if groups.nunique() < 5:
        raise ValueError("At least 5 independent groups are required for this starter's holdout. This is not a clinical sample-size rule.")
    excluded = 0
    if cfg.split_method == "group_random":
        splitter = GroupShuffleSplit(n_splits=1, test_size=cfg.test_size, random_state=cfg.random_state)
        train_idx, test_idx = next(splitter.split(frame, y, groups))
    else:
        if frame[cfg.holdout_column].isna().any():
            raise ValueError("The geographic holdout column must be known for every included row.")
        observed = set(frame[cfg.holdout_column].astype(str))
        if not set(cfg.holdout_values).issubset(observed):
            raise ValueError("A requested held-out location does not occur in the selected cohort.")
        held = frame[cfg.holdout_column].astype(str).isin(cfg.holdout_values)
        held_groups = set(groups[held])
        train_mask = ~held & ~groups.isin(held_groups)
        train_idx, test_idx = np.flatnonzero(train_mask), np.flatnonzero(held)
        excluded = int((~held & groups.isin(held_groups)).sum())
    if not len(train_idx) or not len(test_idx):
        raise ValueError("The split leaves an empty training or test partition.")
    if set(groups.iloc[train_idx]).intersection(set(groups.iloc[test_idx])):
        raise RuntimeError("Internal error: group overlap in split.")
    for name, idx in [("training", train_idx), ("test", test_idx)]:
        counts = y.iloc[idx].value_counts()
        if len(counts) < 2:
            raise ValueError(f"The fixed {name} split contains only one class. Use an adequately represented cohort; do not cherry-pick a seed for a better score.")
    return train_idx, test_idx, excluded


def train_models(frame: pd.DataFrame, cfg: Config, out_dir: str | Path) -> dict[str, Any]:
    cfg.validate()
    _require(frame, [cfg.sample_column, cfg.group_column, cfg.target_column] + cfg.geographic_columns)
    cohort = _scoped(frame, cfg, training=True)
    if len(cohort) < 30:
        raise ValueError("At least 30 labelled rows are required by this software. This minimum is NOT evidence of sufficient training data.")
    for col in {cfg.sample_column, cfg.group_column}:
        if cohort[col].isna().any():
            raise ValueError(f"{col}: every included row needs a de-identified identifier.")
    if cohort[cfg.sample_column].duplicated().any():
        raise ValueError("Repeated isolate IDs for this organism-antibiotic pair. Resolve duplicate/conflicting records before training.")
    y = encode_target(cohort[cfg.target_column])
    if y.nunique() != 2:
        raise ValueError("Both resistant and susceptible examples are needed.")
    train_idx, test_idx, excluded = split_rows(cohort, cfg, y)
    notes = [NOTICE,
             "One pre-specified holdout only; not cross-validation, an external clinical study, or proof of generalization.",
             "The threshold was fixed before evaluation. Do not tune it or select features using this test set.",
             "Global feature importance is model dependence, not a causal mechanism or an explanation of an individual prediction."]
    if min(y.iloc[test_idx].value_counts()) < 20:
        notes.append("Fewer than 20 test examples of at least one class: metrics may be very unstable.")
    if cfg.group_column == cfg.sample_column:
        notes.append("Splitting by isolate only. If one patient has several isolates, use a consistent patient/group identifier instead.")
    if excluded:
        notes.append(f"Excluded {excluded} non-held-location rows sharing groups with the geographic test cohort.")
    # Validate every feature block before writing any model files.
    variants = ["genomic_only", "full"] if cfg.compare_full else ["genomic_only"]
    design_matrices = {}
    for variant in variants:
        X = prepare_features(cohort, cfg, variant)
        empty = X.columns[X.iloc[train_idx].isna().all()].tolist()
        if empty:
            raise ValueError("Features entirely missing in the training partition: " + ", ".join(empty)
                             + ". Remove unsupported features; do not fabricate observations.")
        design_matrices[variant] = X
    out = Path(out_dir)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory is not empty. Use a new directory to preserve previous runs.")
    out.mkdir(parents=True, exist_ok=True)
    report: dict[str, Any] = {
        "format_version": FORMAT_VERSION, "notice": NOTICE,
        "created_utc": datetime.now(timezone.utc).isoformat(), "versions": versions(),
        "config": asdict(cfg), "input_rows": len(frame), "selected_pair_rows": len(cohort),
        "other_pair_rows_not_used": len(frame) - len(cohort), "train_n": len(train_idx), "test_n": len(test_idx),
        "train_groups": int(cohort.iloc[train_idx][cfg.group_column].nunique()),
        "test_groups": int(cohort.iloc[test_idx][cfg.group_column].nunique()),
        "train_class_counts": {str(k): int(v) for k, v in y.iloc[train_idx].value_counts().items()},
        "excluded_shared_holdout_group_rows": excluded, "notes": notes, "models": {},
        "cohort_sha256": hashlib.sha256(cohort.to_csv(index=False).encode()).hexdigest()}
    assignments = cohort[[cfg.sample_column, cfg.group_column]].copy()
    assignments = assignments.loc[:, ~assignments.columns.duplicated()]
    assignments["split"] = "excluded_shared_holdout_group"
    assignments.loc[train_idx, "split"] = "train"
    assignments.loc[test_idx, "split"] = "test"
    assignments.to_csv(out / "split_assignments.csv", index=False)
    heldout = cohort.iloc[test_idx][[cfg.sample_column]].copy().reset_index(drop=True)
    heldout["actual_resistant"] = y.iloc[test_idx].to_numpy()
    slices = []
    for variant in variants:
        X = design_matrices[variant]
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        pipeline = make_pipeline(cfg, variant)
        pipeline.fit(X_train, y.iloc[train_idx])
        scores = _positive_scores(pipeline, X_test)
        heldout[f"{variant}_score_uncalibrated"] = scores
        heldout[f"{variant}_prediction"] = (scores >= cfg.threshold).astype(int)
        details = {"metrics": metrics(y.iloc[test_idx], scores, cfg.threshold),
                   "feature_count": X.shape[1], "input_schema": _input_schema(X_train, cfg, variant)}
        report["models"][variant] = details
        bundle = {"format_version": FORMAT_VERSION, "notice": NOTICE, "versions": versions(),
                  "config": asdict(cfg), "variant": variant, "pipeline": pipeline,
                  "input_schema": details["input_schema"]}
        joblib.dump(bundle, out / f"{variant}.joblib", compress=3)
        if cfg.permutation_repeats:
            importance = permutation_importance(pipeline, X_test, y.iloc[test_idx],
                n_repeats=cfg.permutation_repeats, scoring="roc_auc", random_state=cfg.random_state, n_jobs=1)
            pd.DataFrame({"feature": X.columns, "mean_auc_drop": importance.importances_mean,
                          "std_auc_drop": importance.importances_std}).sort_values("mean_auc_drop", ascending=False).to_csv(
                              out / f"{variant}_global_importance.csv", index=False)
        for col in cfg.geographic_columns:
            locations = cohort.iloc[test_idx][col].fillna("<missing>").astype(str).to_numpy()
            for location in sorted(set(locations)):
                mask = locations == location
                slices.append({"model": variant, "geographic_column": col, "location": location,
                               **metrics(y.iloc[test_idx].to_numpy()[mask], scores[mask], cfg.threshold)})
    if cfg.compare_full:
        report["full_minus_genomic_auc"] = (report["models"]["full"]["metrics"]["roc_auc"]
                                            - report["models"]["genomic_only"]["metrics"]["roc_auc"])
        report["notes"].append("The observed AUC difference is descriptive; it has no confidence interval and may be zero or negative.")
    heldout.to_csv(out / "heldout_predictions.csv", index=False)
    if slices:
        pd.DataFrame(slices).to_csv(out / "geographic_metrics.csv", index=False)
    (out / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    (out / "config.json").write_text(json.dumps(asdict(cfg), indent=2), encoding="utf-8")
    return report


def load_bundle(path: str | Path, *, trusted: bool = False) -> dict[str, Any]:
    """Joblib can execute code while loading. Never load untrusted model files."""
    if not trusted:
        raise ValueError("Only load models you trained locally or obtained from a trusted source. Explicit trust confirmation is required.")
    bundle = joblib.load(path)
    if not isinstance(bundle, dict) or bundle.get("format_version") != FORMAT_VERSION:
        raise ValueError("This is not a supported model bundle.")
    current = versions()
    for package in ("scikit-learn", "numpy", "pandas", "joblib"):
        if bundle["versions"].get(package) != current[package]:
            raise ValueError(f"Version mismatch for {package}. Reuse the training environment or retrain from the original data.")
    return bundle


def predict_rows(bundle: dict[str, Any], frame: pd.DataFrame) -> pd.DataFrame:
    cfg = Config(**bundle["config"])
    cfg.validate()
    cohort = _scoped(frame, cfg, training=False)
    X = prepare_features(cohort, cfg, bundle["variant"])
    flags: list[list[str]] = [[] for _ in range(len(X))]
    for spec in bundle["input_schema"]:
        col = spec["name"]
        missing = X[col].isna()
        for i in np.flatnonzero(missing):
            flags[i].append(f"{col}: missing; training-fitted imputation used")
        if spec["kind"] == "categorical":
            unknown = ~missing & ~X[col].isin(spec["values"])
            for i in np.flatnonzero(unknown):
                flags[i].append(f"{col}: unseen category; outside training categories")
        else:
            outside = ~missing & ((X[col] < spec["min"]) | (X[col] > spec["max"]))
            for i in np.flatnonzero(outside):
                flags[i].append(f"{col}: outside observed training range")
    scores = _positive_scores(bundle["pipeline"], X)
    result = pd.DataFrame({"row_number": np.arange(1, len(X) + 1)})
    if cfg.sample_column in cohort:
        result[cfg.sample_column] = cohort[cfg.sample_column].to_numpy()
    result["organism"] = cfg.organism
    result["antibiotic"] = cfg.antibiotic
    result["model"] = bundle["variant"]
    result["resistance_score_uncalibrated"] = scores
    result["threshold"] = cfg.threshold
    result["prediction"] = np.where(scores >= cfg.threshold, "Predicted resistant", "Predicted susceptible")
    result["missing_feature_count"] = X.isna().sum(axis=1).to_numpy()
    result["input_warnings"] = ["; ".join(x) for x in flags]
    result["research_only"] = True
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    train = commands.add_parser("train", help="Train and evaluate on a fixed group-disjoint holdout")
    train.add_argument("--data", required=True)
    train.add_argument("--config", required=True)
    train.add_argument("--out", required=True)
    predict = commands.add_parser("predict", help="Predict from a locally trained model")
    predict.add_argument("--model", required=True)
    predict.add_argument("--data", required=True)
    predict.add_argument("--out", required=True)
    predict.add_argument("--trust-local-model", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "train":
            cfg = Config(**json.loads(Path(args.config).read_text(encoding="utf-8")))
            report = train_models(read_csv(args.data), cfg, args.out)
            print(NOTICE)
            for name, result in report["models"].items():
                print(name, json.dumps(result["metrics"], indent=2))
            print(f"Saved models and reports to {args.out}")
        else:
            out = Path(args.out)
            if out.exists():
                raise ValueError("Prediction output already exists. Choose a new filename.")
            bundle = load_bundle(args.model, trusted=args.trust_local_model)
            prediction = predict_rows(bundle, read_csv(args.data))
            out.parent.mkdir(parents=True, exist_ok=True)
            prediction.to_csv(out, index=False)
            print(NOTICE)
            print(f"Saved {len(prediction)} predictions to {out}")
        return 0
    except (ValueError, TypeError, OSError, KeyError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
