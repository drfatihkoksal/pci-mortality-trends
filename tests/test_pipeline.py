from __future__ import annotations

import csv
import json
import math
import re
import unittest
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
PRIVATE_COHORT = ROOT / "data" / "private" / "analysis_cohort.csv"
RESULTS = ROOT / "outputs" / "results"
TABLES = ROOT / "outputs" / "tables"
AUDIT = ROOT / "outputs" / "audit"


def keys() -> dict[str, str]:
    with (RESULTS / "key_results.csv").open(encoding="utf-8") as handle:
        return {row["key"]: row["value"] for row in csv.DictReader(handle)}


class PipelineValidation(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cohort = pd.read_csv(PRIVATE_COHORT)
        cls.key = keys()
        cls.audit = json.loads(
            (AUDIT / "data_quality_summary.json").read_text(encoding="utf-8")
        )

    def test_direct_identifiers_are_not_persisted(self) -> None:
        prohibited = {
            "TcKimlik",
            "patientID",
            "HastaAdSoyad",
            "KataterNo",
            "DoktorAdi1",
            "DoktorAdi2",
        }
        self.assertFalse(prohibited.intersection(self.cohort.columns))
        self.assertEqual(self.cohort["study_id"].nunique(), len(self.cohort))
        self.assertTrue(
            self.cohort["study_id"].str.fullmatch(r"P\d{6}").all()
        )

    def test_raw_workbook_is_git_ignored(self) -> None:
        ignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
        self.assertIn("*.xlsx", ignore)

    def test_source_hash_and_counts_match(self) -> None:
        self.assertRegex(self.audit["source_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(int(self.key["source_records"]), self.audit["source_rows"])
        self.assertEqual(int(self.key["unique_patients"]), len(self.cohort))
        self.assertEqual(
            int(self.key["invalid_death_before_index"]),
            int(self.cohort["death_before_index"].sum()),
        )

    def test_fixed_horizon_cohorts_are_mature(self) -> None:
        expected = {
            "primary": (2015, 2024, 365, "primary_n", "primary_events"),
            "early": (2015, 2025, 30, "early_n", "early_events"),
            "long": (2015, 2020, 1825, "long_n", "long_events"),
        }
        valid = self.cohort[
            (self.cohort["death_before_index"] == 0)
            & self.cohort["followup_days"].notna()
            & (self.cohort["followup_days"] >= 0)
        ].copy()
        for label, (first, last, horizon, n_key, event_key) in expected.items():
            with self.subTest(label=label):
                frame = valid[valid["index_year"].between(first, last)].copy()
                event = (frame["dead"] == 1) & (frame["followup_days"] <= horizon)
                mature = (frame["followup_days"] >= horizon) | event
                frame = frame[mature]
                events = ((frame["dead"] == 1) & (frame["followup_days"] <= horizon)).sum()
                self.assertEqual(len(frame), int(self.key[n_key]))
                self.assertEqual(int(events), int(self.key[event_key]))

    def test_primary_flow_reconciles(self) -> None:
        self.assertEqual(
            int(self.key["unique_patients"]),
            int(self.key["outside_primary_period"])
            + int(self.key["invalid_death_before_index_primary"])
            + int(self.key["primary_n"]),
        )

    def test_selected_model_estimates_match_annual_table(self) -> None:
        annual = pd.read_csv(TABLES / "annual_adjusted_mortality.csv")
        selected = annual[
            (annual["outcome"] == "1-year all-cause mortality")
            & (annual["model"] == self.key["selected_model"])
        ].set_index("index_year")
        self.assertAlmostEqual(
            100 * selected.loc[2015, "risk"],
            float(self.key["primary_start_adjusted_percent"]),
            places=9,
        )
        self.assertAlmostEqual(
            100 * selected.loc[2024, "risk"],
            float(self.key["primary_end_adjusted_percent"]),
            places=9,
        )
        self.assertTrue(selected["risk"].between(0, 1).all())
        self.assertTrue((selected["lower"] <= selected["risk"]).all())
        self.assertTrue((selected["risk"] <= selected["upper"]).all())

    def test_manuscript_has_no_editorial_notes(self) -> None:
        paper = (ROOT / "manuscript" / "paper.qmd").read_text(encoding="utf-8")
        supplement = (ROOT / "manuscript" / "supplement.md").read_text(encoding="utf-8")
        for text in (paper, supplement):
            self.assertNotRegex(text, r"(?i)internal working draft|before (journal )?submission|working analysis|required\]")
        self.assertIn("[@dawson2021]", paper)
        self.assertRegex(paper, r"22,171 patients")

    def test_bibliography_has_verified_dois(self) -> None:
        bibliography = (ROOT / "manuscript" / "references.bib").read_text(
            encoding="utf-8"
        )
        entries = re.findall(r"@article\{.*?\n\}", bibliography, flags=re.S)
        self.assertGreaterEqual(len(entries), 10)
        for entry in entries:
            with self.subTest(entry=entry[:40]):
                self.assertRegex(entry, r"doi\s*=\s*\{10\.\S+?\}")

    def test_procedure_days_are_deidentified(self) -> None:
        days = pd.read_csv(ROOT / "data" / "private" / "procedure_days.csv")
        prohibited = {"TcKimlik", "patientID", "HastaAdSoyad", "KataterNo", "DoktorAdi1", "DoktorAdi2"}
        self.assertFalse(prohibited.intersection(days.columns))
        self.assertEqual(len(days), self.audit["procedure_days"])
        self.assertTrue(set(days["study_id"]).issuperset(set(self.cohort["study_id"])))

    def test_sensitivity_primary_row_matches_key_results(self) -> None:
        sens = pd.read_csv(TABLES / "sensitivity_analyses.csv")
        first = sens[sens["outcome"] == "1-year all-cause mortality"].iloc[0]
        self.assertTrue(first["analysis"].startswith("Primary"))
        self.assertAlmostEqual(
            100 * first["start_risk"], float(self.key["primary_start_adjusted_percent"]), places=9
        )
        self.assertAlmostEqual(
            100 * first["rd_end_start"], float(self.key["primary_risk_difference_points"]), places=9
        )
        self.assertTrue((sens["rd_peak_start_lower"] <= sens["rd_peak_start"]).all())

    def test_every_citation_is_in_bibliography(self) -> None:
        paper = (ROOT / "manuscript" / "paper.qmd").read_text(encoding="utf-8")
        bibliography = (ROOT / "manuscript" / "references.bib").read_text(encoding="utf-8")
        cited = set(re.findall(r"@([a-z]+\d{4})", paper))
        keys = set(re.findall(r"@article\{([^,]+),", bibliography))
        self.assertTrue(cited)
        self.assertFalse(cited - keys, f"missing: {sorted(cited - keys)}")

    def test_required_outputs_exist_and_are_nonempty(self) -> None:
        paths = [
            TABLES / "annual_crude_mortality.csv",
            TABLES / "annual_adjusted_mortality.csv",
            TABLES / "model_contrasts.csv",
            TABLES / "sensitivity_analyses.csv",
            TABLES / "source_bridge_2022.csv",
            TABLES / "period_mortality_2015_2018_cohorts.csv",
            ROOT / "outputs" / "figures" / "figureS1_sensitivity_trends.png",
            ROOT / "outputs" / "figures" / "figure1_cohort_flow.png",
            ROOT / "outputs" / "figures" / "figure3_adjusted_mortality_trends.png",
            ROOT / "manuscript" / "paper.qmd",
            ROOT / "manuscript" / "supplement.md",
            ROOT / "submission" / "cover_letter.md",
        ]
        for path in paths:
            with self.subTest(path=path):
                self.assertTrue(path.exists())
                self.assertGreater(path.stat().st_size, 100)


if __name__ == "__main__":
    unittest.main()
