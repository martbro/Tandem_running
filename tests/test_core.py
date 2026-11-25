"""Regression cases for the shared classification path."""

import unittest

from acmg_snv.criteria import score_class
from acmg_snv.engine import Classifier, load_ps3
from acmg_snv.model import Variant


SNV = Variant.parse("chr1:123:A:G")
BASE = {
    "gene": "GENE_DEMO", "transcript": "NM_TEST.1", "mane": "MANE Select",
    "consequence": "missense_variant", "source": "annotation simulée",
}
CONTEXT = {"inheritance": "AD", "disease": "Maladie simulée"}


def override(code, strength, source):
    return {
        "code": code, "strength": strength, "reason": "Preuve simulée",
        "source": source, "evidence_ids": [source],
    }


def row(report, code):
    return next(item for item in report["criteria"] if item["code"] == code)


class CoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = Classifier()
        cls.addClassCleanup(cls.engine.close)

    def check(self, annotation=None, overrides=None, **context):
        return self.engine.classify(
            SNV, CONTEXT | context, BASE | (annotation or {}), overrides,
        )

    def test_2026_score_edges(self):
        values = [-8, -7, -6, -1, 0, 5, 6, 9, 10]
        expected = ["B", "B", "LB", "LB", "VUS", "VUS", "LP", "LP", "P"]
        self.assertEqual([score_class(value) for value in values], expected)

    def test_single_criterion_requires_corroboration(self):
        data = [override("PVS1", "very_strong", "preuve-1")]
        report = self.check(overrides=data)
        self.assertEqual(report["score"], 8)
        self.assertEqual(report["classification"], "VUS")

    def test_bs1_standalone(self):
        data = [override("BS1", "strong", "preuve-2")]
        standard = self.check(overrides=data)
        standalone = self.check(overrides=data, high_penetrance=True)
        self.assertEqual(standard["classification"], "VUS")
        self.assertEqual(standalone["classification"], "LB")

    def test_ba1_standalone(self):
        data = [override("BA1", "standalone", "preuve-2")]
        report = self.check(overrides=data)
        self.assertEqual(report["classification"], "B")
        self.assertIsNone(row(report, "BA1")["points"])

    def test_negative_six_is_likely_benign(self):
        report = self.check(overrides=[
            override("BS2", "strong", "témoins"),
            override("BS3", "moderate", "essais"),
        ])
        self.assertEqual(report["classification"], "LB")
        self.assertEqual(report["score"], -6)

    def test_conflict_requires_expert(self):
        evidence = [
            override("PS3", "strong", "essai"),
            override("BS2", "strong", "population"),
        ]
        report = self.check(overrides=evidence)
        self.assertEqual(report["status"], "review_required")
        self.assertIsNone(report["classification"])
        accepted = self.check(
            overrides=evidence, reviewed=True,
            conflict_resolution="Les observations indépendantes sont retenues.",
        )
        self.assertEqual(accepted["classification"], "VUS")

    def test_reused_evidence_is_quarantined(self):
        report = self.check(overrides=[
            override("PS3", "strong", "identique"),
            override("PS4", "strong", "identique"),
        ])
        self.assertEqual(report["status"], "review_required")
        self.assertIsNone(report["classification"])
        self.assertTrue(all(
            row(report, code)["status"] != "applied" for code in ("PS3", "PS4")
        ))

    def test_functional_file_takes_priority(self):
        text = "variant\tstrength\treason\nchr1:123:A:G\tstrong\tEssai.\n"
        annotation = BASE | {"functional_ps3": {
            "validated": True, "strength": "supporting", "reason": "papier",
            "source": "PMID simulé", "evidence_id": "papier",
        }}
        report = self.engine.classify(
            SNV, CONTEXT, annotation, ps3=load_ps3(text=text),
        )
        self.assertEqual(row(report, "PS3")["points"], 4)
        self.assertEqual(row(report, "PS3")["origin"], "laboratory_ps3")

    def test_non_scoped_and_unknown_context(self):
        outside = self.check(annotation={"consequence": "intergenic_variant"})
        absent = self.engine.classify(SNV, {})
        self.assertEqual(outside["status"], "out_of_scope")
        self.assertEqual(absent["status"], "needs_context")

    def test_vcep_prevents_generic_scoring(self):
        report = self.check(
            annotation={"gene": "MYH7", "revel": 0.97},
            disease="MONDO:0005045",
        )
        self.assertEqual(report["status"], "review_required")
        self.assertEqual(report["vcep"]["matches"][0]["id"], "GN002")
        self.assertEqual(row(report, "PP3")["status"], "not_evaluable")

    def test_denovo_criteria_must_be_aggregated(self):
        report = self.check(overrides=[
            override("PS2", "strong", "observation A"),
            override("PM6", "strong", "observation B"),
        ])
        self.assertEqual(report["status"], "review_required")
        self.assertNotEqual(row(report, "PM6")["status"], "applied")

    def test_vcep_incompatible_inheritance_does_not_fall_back(self):
        report = self.check(
            annotation={"gene": "MYH7", "revel": 0.97},
            disease="MONDO:0005045", inheritance="AR",
            vcep_exclusion_reason="Exclusion simulée",
        )
        self.assertEqual(report["status"], "needs_context")

    def test_vcep_manual_confirmation_needs_evidence(self):
        report = self.check(
            annotation={"gene": "MYH7"}, disease="MONDO:0005045",
            vcep_reviewed=True, reviewed=True,
        )
        self.assertEqual(report["status"], "not_evaluable")

    def test_row_without_gene_enriches_transcript(self):
        annotation = [BASE, {
            "population_ac": 0, "population_an": 80_000,
            "population_group": "groupe simulé", "population_source": "gnomAD",
            "population_pass": True, "population_eligible": True,
            "population_absence_verified": True,
        }]
        report = self.engine.classify(SNV, CONTEXT, annotation)
        self.assertEqual(row(report, "PM2")["points"], 2)

    def test_manual_changes_require_a_reason(self):
        with self.assertRaisesRegex(ValueError, "Justification"):
            self.check(overrides=[{"code": "PP4", "strength": "supporting"}])


if __name__ == "__main__":
    unittest.main()
