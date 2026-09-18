"""Software tests for the data-preparation, metadata, phenotype and prediction tooling.

These check behaviour the rest of the project depends on being correct: that
feature column names are reproducible across processes, that locus annotation
includes promoter margins, that ambiguous laboratory results are dropped rather
than recoded, that free-text metadata placeholders are not mistaken for data,
and that the standalone predictor reproduces the training pipeline exactly.

They use small synthetic fixtures with deliberately non-biological values where
the biology does not matter. They establish software correctness, not clinical
validity or predictive performance.
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import build_phenotype_table as bpt
import fetch_sample_metadata as fsm
import predict_resistance as pr
import prepare_variants as pv

REPO = Path(__file__).resolve().parent.parent


class FeatureNamingTests(unittest.TestCase):
    def test_short_allele_kept_verbatim(self):
        self.assertEqual(pv._shorten("ACGT"), "ACGT")

    def test_long_allele_is_shortened_with_a_digest(self):
        name = pv._shorten("ACGTACGTACGTACGTACGT")
        self.assertTrue(name.startswith("ACGTACGTACGT"))
        self.assertIn("x", name)

    def test_digest_is_stable_across_processes(self):
        """Built-in hash() is salted per process; the digest must not be."""
        code = ("import sys; sys.path.insert(0, %r); import prepare_variants as p; "
                "print(p._shorten('ACGTACGTACGTACGTACGT'))" % str(REPO))
        first = subprocess.run([sys.executable, "-c", code], capture_output=True,
                               text=True, check=True).stdout.strip()
        second = subprocess.run([sys.executable, "-c", code], capture_output=True,
                                text=True, check=True).stdout.strip()
        self.assertEqual(first, second)
        self.assertEqual(first, pv._shorten("ACGTACGTACGTACGTACGT"))

    def test_distinct_long_alleles_get_distinct_names(self):
        a = pv._shorten("ACGTACGTACGTACGTACGTAA")
        b = pv._shorten("ACGTACGTACGTACGTACGTGG")
        self.assertNotEqual(a, b)

    def test_column_name_includes_gene_and_position(self):
        self.assertEqual(pv.column_name("rpoB", "761155_C_T"), "g_rpoB_761155_C_T")

    def test_unannotated_position_is_labelled_intergenic(self):
        self.assertTrue(pv.column_name("", "1234_A_G").startswith("g_intergenic_"))


class LocusAnnotationTests(unittest.TestCase):
    def setUp(self):
        self.index = pv.LocusIndex(
            [{"gene": "geneA", "start": 1000, "end": 2000, "drugs": "drugX", "note": ""},
             {"gene": "geneB", "start": 5000, "end": 6000, "drugs": "drugY", "note": ""}],
            margin=100)

    def test_position_inside_a_locus_is_annotated(self):
        self.assertEqual(list(self.index.gene_of(np.array([1500]))), ["geneA"])

    def test_margin_captures_upstream_promoter_positions(self):
        # 950 is outside geneA's coding range but inside the 100-base margin.
        self.assertEqual(list(self.index.gene_of(np.array([950]))), ["geneA"])

    def test_position_outside_every_locus_is_blank(self):
        self.assertEqual(list(self.index.gene_of(np.array([3000]))), [""])

    def test_first_matching_locus_wins(self):
        overlapping = pv.LocusIndex(
            [{"gene": "first", "start": 10, "end": 100, "drugs": "d", "note": ""},
             {"gene": "second", "start": 50, "end": 200, "drugs": "d", "note": ""}],
            margin=0)
        self.assertEqual(list(overlapping.gene_of(np.array([60]))), ["first"])


class PhenotypeLoadingTests(unittest.TestCase):
    def _write(self, text: str) -> Path:
        path = Path(self.tmp.name) / "pheno.csv"
        path.write_text(text, encoding="utf-8")
        return path

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def test_r_and_s_become_one_and_zero(self):
        path = self._write("SAMPLE,DRUG\nA,R\nB,S\n")
        frame = pv.load_phenotypes(path, "SAMPLE", "DRUG")
        self.assertEqual(dict(zip(frame.isolate_id, frame.resistant)), {"A": 1, "B": 0})

    def test_intermediate_is_dropped_not_recoded(self):
        path = self._write("SAMPLE,DRUG\nA,R\nB,I\nC,S\n")
        frame = pv.load_phenotypes(path, "SAMPLE", "DRUG")
        self.assertEqual(set(frame.isolate_id), {"A", "C"})
        self.assertEqual(frame.attrs["dropped_unusable_labels"], 1)

    def test_blank_and_unknown_are_dropped(self):
        path = self._write("SAMPLE,DRUG\nA,R\nB,\nC,unknown\n")
        frame = pv.load_phenotypes(path, "SAMPLE", "DRUG")
        self.assertEqual(list(frame.isolate_id), ["A"])

    def test_conflicting_duplicate_isolate_aborts(self):
        path = self._write("SAMPLE,DRUG\nA,R\nA,S\n")
        with self.assertRaises(SystemExit):
            pv.load_phenotypes(path, "SAMPLE", "DRUG")

    def test_missing_drug_column_aborts(self):
        path = self._write("SAMPLE,DRUG\nA,R\n")
        with self.assertRaises(SystemExit):
            pv.load_phenotypes(path, "SAMPLE", "NOT_A_COLUMN")


class MetadataNormalizationTests(unittest.TestCase):
    def setUp(self):
        rules = json.loads((REPO / "resources" / "metadata_normalization.json")
                           .read_text(encoding="utf-8"))
        self.normalizer = fsm.Normalizer(rules)

    def test_literal_missing_is_treated_as_absent(self):
        self.assertTrue(self.normalizer.is_null("missing"))
        self.assertTrue(self.normalizer.is_null("Not applicable"))
        self.assertTrue(self.normalizer.is_null(""))

    def test_a_real_country_is_not_null(self):
        self.assertFalse(self.normalizer.is_null("South Africa"))

    def test_facility_names_are_recognised_and_not_countries(self):
        self.assertTrue(self.normalizer.is_facility("Catherine Booth"))
        self.assertFalse(self.normalizer.is_facility("Peru"))

    def test_specimen_grouping(self):
        self.assertEqual(self.normalizer.specimen_category("sputum"), "respiratory")
        self.assertEqual(self.normalizer.specimen_category("cerebrospinal fluid"),
                         "extrapulmonary")

    def test_non_anatomical_source_is_unspecified_not_respiratory(self):
        self.assertEqual(self.normalizer.specimen_category("patient"), "unspecified")
        self.assertEqual(self.normalizer.specimen_category("culture"), "unspecified")

    def test_null_specimen_stays_blank(self):
        self.assertEqual(self.normalizer.specimen_category("missing"), "")

    def test_country_and_region_split(self):
        self.assertEqual(fsm.split_country("United Kingdom: Midlands"),
                         ("United Kingdom", "Midlands"))
        self.assertEqual(fsm.split_country("Peru"), ("Peru", ""))
        self.assertEqual(fsm.split_country(""), ("", ""))

    def test_collection_year_parsing_rejects_nonsense(self):
        self.assertEqual(fsm.collection_year("2004-11-01"), "2004")
        self.assertEqual(fsm.collection_year("2004"), "2004")
        self.assertEqual(fsm.collection_year("not collected"), "")
        self.assertEqual(fsm.collection_year("1780"), "")


class PhenotypeMergeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def test_ncbi_ast_parsing(self):
        path = self.dir / "ncbi.tsv"
        path.write_text(
            "Run\tAST_phenotypes\n"
            "SRR1\trifampin=R,isoniazid=S,ethambutol=I\n"
            "SRR2,SRR3\tpyrazinamide=S\n"
            "SRR4\tNULL\n", encoding="utf-8")
        frame = bpt.load_ncbi(path)
        pairs = {(r.isolate_id, r.drug): r.label for r in frame.itertuples()}
        self.assertEqual(pairs[("SRR1", "RIFAMPICIN")], "R")
        self.assertEqual(pairs[("SRR1", "ISONIAZID")], "S")
        # Intermediate is not a class being modelled, so it is dropped.
        self.assertNotIn(("SRR1", "ETHAMBUTOL"), pairs)
        # One row naming several runs applies to each of them.
        self.assertEqual(pairs[("SRR2", "PYRAZINAMIDE")], "S")
        self.assertEqual(pairs[("SRR3", "PYRAZINAMIDE")], "S")
        self.assertNotIn(("SRR4", "PYRAZINAMIDE"), pairs)

    def test_cryptic_quality_filter(self):
        path = self.dir / "cryptic.csv"
        path.write_text(
            "ENA_RUN,RIF_BINARY_PHENOTYPE,RIF_PHENOTYPE_QUALITY\n"
            "ERR1,R,HIGH\nERR2,S,LOW\nERR3,R,MEDIUM\n", encoding="utf-8")
        strict = bpt.load_cryptic(path, {"HIGH", "MEDIUM"})
        self.assertEqual(set(strict.isolate_id), {"ERR1", "ERR3"})
        permissive = bpt.load_cryptic(path, set())
        self.assertEqual(set(permissive.isolate_id), {"ERR1", "ERR2", "ERR3"})


class PredictorTests(unittest.TestCase):
    """The standalone predictor must rebuild training columns exactly."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def test_column_to_key_uses_the_dictionary(self):
        dictionary = pd.DataFrame({
            "feature": ["g_rpoB_761155_C_T", "g_katG_2155168_C_G"],
            "position": ["761155", "2155168"],
            "ref": ["C", "C"], "alt": ["T", "G"]})
        mapping = pr.column_to_key(["g_rpoB_761155_C_T"], dictionary)
        self.assertEqual(mapping["g_rpoB_761155_C_T"], "761155_C_T")

    def test_unknown_feature_aborts_rather_than_scoring_blind(self):
        dictionary = pd.DataFrame({"feature": ["g_a_1_C_T"], "position": ["1"],
                                   "ref": ["C"], "alt": ["T"]})
        with self.assertRaises(SystemExit):
            pr.column_to_key(["g_missing_9_A_G"], dictionary)

    def test_matrix_switches_on_only_known_variants(self):
        variants = self.dir / "v.csv"
        variants.write_text(
            "SAMPLE,CHROM,POS,REF,ALT\n"
            "S1,C,761155,C,T\n"      # known
            "S1,C,999999,A,G\n"      # unknown, must be ignored
            "S2,C,111111,T,A\n",     # unknown only
            encoding="utf-8")
        key_of_column = {"g_rpoB_761155_C_T": "761155_C_T",
                         "g_katG_2155168_C_G": "2155168_C_G"}
        matrix, calls = pr.build_matrix(variants, key_of_column, progress=False)
        self.assertEqual(list(matrix["isolate_id"]), ["S1", "S2"])
        self.assertEqual(matrix.loc[0, "g_rpoB_761155_C_T"], 1)
        self.assertEqual(matrix.loc[0, "g_katG_2155168_C_G"], 0)
        # Every supplied call is counted, including ones with no column.
        self.assertEqual(list(calls), [2, 1])
        # An isolate with no recognised feature is all zeros, and callers must
        # flag it rather than read it as susceptible.
        self.assertEqual(int(matrix.loc[1, list(key_of_column)].sum()), 0)

    def test_missing_input_columns_abort(self):
        variants = self.dir / "bad.csv"
        variants.write_text("SAMPLE,POS\nS1,1\n", encoding="utf-8")
        with self.assertRaises(SystemExit):
            pr.build_matrix(variants, {"g_a_1_C_T": "1_C_T"}, progress=False)

    def test_predict_requires_explicit_model_trust(self):
        with self.assertRaises(SystemExit):
            pr.parse_args(["--variants", "v.csv", "--out", "o.csv"])

    def test_threshold_must_be_a_probability(self):
        with self.assertRaises(SystemExit):
            pr.parse_args(["--variants", "v.csv", "--out", "o.csv",
                           "--trust-local-models", "--threshold", "1.5"])


if __name__ == "__main__":
    unittest.main()
