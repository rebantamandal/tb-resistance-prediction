"""Software tests only. Generated observations are NOT biological training data."""
import io
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from amr import (Config, encode_target, load_bundle, metrics, predict_rows,
                 prepare_features, read_csv, split_rows, train_models)


def dummy_frame(n=240):
    """Abstract synthetic fixture, deliberately not a real organism or drug."""
    rng = np.random.default_rng(170)
    return pd.DataFrame({
        "isolate_id": [f"DUMMY_I_{i:04d}" for i in range(n)],
        "patient_id": [f"DUMMY_P_{i // 2:04d}" for i in range(n)],
        "organism": ["DEMO_NONBIOLOGICAL"] * n,
        "antibiotic": ["DEMO_NOT_A_DRUG"] * n,
        "genomic_feature_1": rng.integers(0, 2, n),
        "genomic_feature_2": rng.integers(0, 2, n),
        "age": rng.integers(20, 80, n),
        "prior_antibiotic_exposure": rng.integers(0, 2, n),
        "hospitalization_days_before_sample": rng.integers(0, 20, n),
        "infection_site": rng.choice(["SITE_A", "SITE_B"], n),
        "country": ["PLACE_A" if i < n // 2 else "PLACE_B" for i in range(n)],
        "city": ["CITY_A" if i < n // 2 else "CITY_B" for i in range(n)],
        "resistant": np.arange(n) % 2,
    })


def dummy_config(**overrides):
    base = Config(organism="DEMO_NONBIOLOGICAL", antibiotic="DEMO_NOT_A_DRUG",
                  genomic_columns=["genomic_feature_1", "genomic_feature_2"],
                  clinical_numeric_columns=["age", "prior_antibiotic_exposure", "hospitalization_days_before_sample"],
                  clinical_categorical_columns=["infection_site"],
                  geographic_columns=["country", "city"], n_estimators=30)
    return replace(base, **overrides)


class ModelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        cls.frame = dummy_frame()
        cls.cfg = dummy_config()
        cls.report = train_models(cls.frame, cls.cfg, cls.root / "main")
        cls.bundle = load_bundle(cls.root / "main" / "full.joblib", trusted=True)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def new_input(self):
        return self.frame.head(3).drop(columns="resistant").copy()

    def test_two_models_same_holdout(self):
        out = pd.read_csv(self.root / "main" / "heldout_predictions.csv")
        self.assertEqual(len(out), self.report["test_n"])
        self.assertTrue(out[["genomic_only_score_uncalibrated", "full_score_uncalibrated"]].notna().all().all())

    def test_group_disjoint(self):
        split = pd.read_csv(self.root / "main" / "split_assignments.csv")
        train = set(split.loc[split.split.eq("train"), "patient_id"])
        test = set(split.loc[split.split.eq("test"), "patient_id"])
        self.assertFalse(train.intersection(test))

    def test_serialization_and_prediction(self):
        result = predict_rows(self.bundle, self.new_input())
        self.assertEqual(len(result), 3)
        self.assertTrue(result.resistance_score_uncalibrated.between(0, 1).all())
        self.assertTrue(result.research_only.all())

    def test_target_not_needed_for_prediction(self):
        result = predict_rows(self.bundle, self.new_input())
        self.assertEqual(len(result), 3)

    def test_missing_feature_rejected(self):
        with self.assertRaisesRegex(ValueError, "Missing required columns"):
            predict_rows(self.bundle, self.new_input().drop(columns="genomic_feature_1"))

    def test_missing_value_flagged(self):
        data = self.new_input().astype(object)
        data.loc[0, "age"] = np.nan
        result = predict_rows(self.bundle, data)
        self.assertIn("age: missing", result.loc[0, "input_warnings"])

    def test_unseen_category_flagged(self):
        data = self.new_input()
        data.loc[0, "country"] = "UNSEEN_PLACE"
        result = predict_rows(self.bundle, data)
        self.assertIn("country: unseen category", result.loc[0, "input_warnings"])

    def test_out_of_range_flagged(self):
        data = self.new_input()
        data.loc[0, "age"] = 120
        result = predict_rows(self.bundle, data)
        self.assertIn("age: outside", result.loc[0, "input_warnings"])

    def test_wrong_organism_rejected(self):
        data = self.new_input()
        data.loc[0, "organism"] = "OTHER"
        with self.assertRaisesRegex(ValueError, "trained organism-antibiotic pair"):
            predict_rows(self.bundle, data)

    def test_wrong_antibiotic_rejected(self):
        data = self.new_input()
        data.loc[0, "antibiotic"] = "OTHER"
        with self.assertRaises(ValueError):
            predict_rows(self.bundle, data)

    def test_intermediate_target_not_reclassified(self):
        with self.assertRaisesRegex(ValueError, "NOT silently reclassified"):
            encode_target(pd.Series(["R", "S", "I"]))

    def test_unknown_target_not_reclassified(self):
        with self.assertRaises(ValueError):
            encode_target(pd.Series([0, 1, -1]))

    def test_missing_target_rejected(self):
        with self.assertRaises(ValueError):
            encode_target(pd.Series([0, 1, np.nan]))

    def test_duplicate_isolate_rejected(self):
        data = self.frame.copy()
        data.loc[1, "isolate_id"] = data.loc[0, "isolate_id"]
        with self.assertRaisesRegex(ValueError, "Repeated isolate"):
            train_models(data, self.cfg, self.root / "duplicate")

    def test_binary_genomic_validation(self):
        data = self.new_input()
        data.loc[0, "genomic_feature_1"] = -1
        with self.assertRaisesRegex(ValueError, "0, 1, or blank"):
            predict_rows(self.bundle, data)

    def test_entirely_missing_column_rejected(self):
        data = self.frame.copy()
        data["age"] = np.nan
        with self.assertRaisesRegex(ValueError, "entirely missing"):
            train_models(data, self.cfg, self.root / "all_missing")

    def test_no_observed_predictors_rejected(self):
        data = self.new_input().astype(object)
        for col in self.cfg.genomic_columns:
            data[col] = np.nan
        with self.assertRaisesRegex(ValueError, "no observed predictors"):
            prepare_features(data, self.cfg, "genomic_only")

    def test_negative_age_rejected(self):
        data = self.new_input()
        data.loc[0, "age"] = -1
        with self.assertRaisesRegex(ValueError, "negative"):
            predict_rows(self.bundle, data)

    def test_preprocessing_fitted_on_training_only(self):
        split = pd.read_csv(self.root / "main" / "split_assignments.csv")
        ids = set(split.loc[split.split.eq("train"), "isolate_id"])
        expected = self.frame.loc[self.frame.isolate_id.isin(ids), "age"].median()
        imputer = self.bundle["pipeline"].named_steps["prepare"].named_transformers_["clinical_numeric"]
        self.assertEqual(imputer.statistics_[0], expected)

    def test_target_cannot_be_feature(self):
        cfg = replace(self.cfg, genomic_columns=["resistant"])
        with self.assertRaisesRegex(ValueError, "cannot be predictor"):
            cfg.validate()

    def test_full_requires_matched_feature_blocks(self):
        cfg = replace(self.cfg, geographic_columns=[])
        with self.assertRaisesRegex(ValueError, "clinical AND geographic"):
            cfg.validate()

    def test_identifier_preservation_in_csv(self):
        data = read_csv(io.StringIO("id,country\n000012,NA\n"))
        self.assertEqual(data.loc[0, "id"], "000012")
        self.assertEqual(data.loc[0, "country"], "NA")

    def test_duplicate_headers_rejected(self):
        with self.assertRaisesRegex(ValueError, "unique"):
            read_csv(io.StringIO("a,a\n1,2\n"))

    def test_geographic_holdout(self):
        cfg = replace(self.cfg, split_method="geography_holdout", holdout_column="country", holdout_values=["PLACE_B"])
        tr, te, excluded = split_rows(self.frame, cfg, encode_target(self.frame.resistant))
        self.assertTrue(self.frame.iloc[te].country.eq("PLACE_B").all())
        self.assertTrue(self.frame.iloc[tr].country.eq("PLACE_A").all())
        self.assertEqual(excluded, 0)

    def test_geographic_shared_group_excluded_from_training(self):
        data = self.frame.copy()
        data.loc[0, "patient_id"] = data.loc[120, "patient_id"]
        cfg = replace(self.cfg, split_method="geography_holdout", holdout_column="country", holdout_values=["PLACE_B"])
        tr, te, excluded = split_rows(data, cfg, encode_target(data.resistant))
        self.assertGreater(excluded, 0)
        self.assertFalse(set(data.iloc[tr].patient_id) & set(data.iloc[te].patient_id))

    def test_single_class_slice_metrics(self):
        result = metrics([0, 0], [0.2, 0.4], 0.5)
        self.assertIsNone(result["roc_auc"])
        self.assertIsNone(result["sensitivity"])
        self.assertEqual(result["specificity"], 1.0)

    def test_loading_untrusted_model_rejected(self):
        with self.assertRaisesRegex(ValueError, "trust"):
            load_bundle(self.root / "main" / "full.joblib")

    def test_refuses_overwriting_run(self):
        with self.assertRaisesRegex(ValueError, "not empty"):
            train_models(self.frame, self.cfg, self.root / "main")

    def test_genomic_only_mode(self):
        cfg = replace(self.cfg, compare_full=False, clinical_numeric_columns=[],
                      clinical_categorical_columns=[], geographic_columns=[])
        report = train_models(self.frame, cfg, self.root / "baseline")
        self.assertEqual(list(report["models"]), ["genomic_only"])

    def test_global_importance_file(self):
        cfg = replace(self.cfg, compare_full=False, permutation_repeats=2)
        train_models(self.frame, cfg, self.root / "importance")
        importance = pd.read_csv(self.root / "importance" / "genomic_only_global_importance.csv")
        self.assertEqual(set(importance.feature), set(self.cfg.genomic_columns))

    def test_report_is_strict_json(self):
        report = json.loads((self.root / "main" / "report.json").read_text())
        self.assertIn("RESEARCH ONLY", report["notice"])
        self.assertEqual(report["config"]["threshold"], 0.5)


if __name__ == "__main__":
    unittest.main()
