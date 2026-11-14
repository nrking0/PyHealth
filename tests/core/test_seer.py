"""Tests for SEERDataset, prepare_seer and parse_sas_layout.

All data are synthetic. The shared fixture in ``test-resources/seer`` holds
seven invented patients (eight tumors) as a fixed-width ASCII file with its
SAS layout, a SEER*Stat style case list, and the expected ``seer.csv`` for
the lung 2004-2009 cohort (five patients, six tumors). Tests copy what they
need into temporary directories.
"""

import gzip
import itertools
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from pyhealth.datasets import SEERDataset
from pyhealth.datasets.seer import (
    LUNG_SITE_RECODE,
    SEER_COLUMNS,
    SEER_FIELDS,
    SEER_GEOGRAPHY_COLUMNS,
    SEER_OUTCOME_COLUMNS,
    parse_sas_layout,
    prepare_seer,
)

FIXTURE = Path(__file__).parent.parent.parent / "test-resources" / "seer"
SAS = FIXTURE / "read.seer.synthetic.sas"
INCIDENCE = FIXTURE / "incidence"
CASES = FIXTURE / "CASES.csv"
COHORT = (2004, 2009)
LUNG = [LUNG_SITE_RECODE]


def read_str(path: Path) -> pd.DataFrame:
    """Reads a CSV with every value as a string and blanks as ''."""
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def labelled_export(directory: Path, drop: tuple[str, ...] = ()) -> Path:
    """Writes the fixture cohort as a SEER*Stat style export (display names)."""
    table = read_str(FIXTURE / "seer.csv").drop(columns=["case_id", *drop])
    labels = {field.name: field.labels[0] for field in SEER_FIELDS}
    path = directory / "export.csv"
    table.rename(columns=labels).to_csv(path, index=False)
    return path


class TempDirTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)


class TestColumns(unittest.TestCase):
    """The canonical columns agree with the YAML config."""

    def test_yaml_attributes_match_columns(self):
        config_path = (
            Path(__file__).parent.parent.parent
            / "pyhealth"
            / "datasets"
            / "configs"
            / "seer.yaml"
        )
        table = yaml.safe_load(config_path.read_text())["tables"]["seer"]
        self.assertEqual(table["patient_id"], "patient_id")
        self.assertEqual(
            table["attributes"], [c for c in SEER_COLUMNS if c != "patient_id"]
        )

    def test_column_groups(self):
        self.assertEqual(SEER_GEOGRAPHY_COLUMNS, ("registry", "state_county"))
        self.assertIn("survival_months", SEER_OUTCOME_COLUMNS)
        self.assertIn("cause_of_death_class", SEER_OUTCOME_COLUMNS)
        self.assertFalse(set(SEER_GEOGRAPHY_COLUMNS) & set(SEER_OUTCOME_COLUMNS))
        self.assertEqual(len(SEER_COLUMNS), len(set(SEER_COLUMNS)))

    def test_fixture_has_canonical_columns(self):
        self.assertEqual(tuple(read_str(FIXTURE / "seer.csv").columns), SEER_COLUMNS)


class TestParseSasLayout(TempDirTestCase):
    def test_fixture_layout(self):
        columns = parse_sas_layout(SAS)
        self.assertEqual(len(columns), 30)
        first = columns[0]
        self.assertEqual(
            (first.name, first.start, first.width, first.label),
            ("PUBCSNUM", 1, 8, "Patient ID"),
        )
        # Positions are contiguous in the synthetic layout.
        for before, after in itertools.pairwise(columns):
            self.assertEqual(after.start, before.start + before.width)
        self.assertIn("RX_SUMM_SURG", [c.name for c in columns])

    def test_duplicate_variable_raises(self):
        path = self.tmp / "dup.sas"
        path.write_text("@ 1 PUBCSNUM $char8.\n@ 9 PUBCSNUM $char8.\n")
        with self.assertRaises(ValueError):
            parse_sas_layout(path)

    def test_no_columns_raises(self):
        path = self.tmp / "empty.sas"
        path.write_text("data seer;\nrun;\n")
        with self.assertRaises(ValueError):
            parse_sas_layout(path)

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            parse_sas_layout(self.tmp / "nope.sas")


class TestPrepareFromAscii(TempDirTestCase):
    def prepare(self, **kwargs) -> pd.DataFrame:
        kwargs.setdefault("incidence_path", INCIDENCE)
        kwargs.setdefault("sas_path", SAS)
        return read_str(prepare_seer(self.tmp / "out", **kwargs))

    def test_year_and_site_filter_matches_fixture(self):
        table = self.prepare(year_range=COHORT, site_recodes=LUNG)
        pd.testing.assert_frame_equal(table, read_str(FIXTURE / "seer.csv"))

    def test_case_list_matches_fixture(self):
        table = self.prepare(cases_path=CASES)
        pd.testing.assert_frame_equal(table, read_str(FIXTURE / "seer.csv"))

    def test_gzipped_file(self):
        directory = self.tmp / "gz"
        directory.mkdir()
        with (
            open(INCIDENCE / "SYNTH.TXT", "rb") as src,
            gzip.open(directory / "SYNTH.TXT.gz", "wb") as dst,
        ):
            dst.write(src.read())
        table = self.prepare(
            incidence_path=directory, year_range=COHORT, site_recodes=LUNG
        )
        pd.testing.assert_frame_equal(table, read_str(FIXTURE / "seer.csv"))

    def test_no_filter_keeps_all_cases(self):
        table = self.prepare()
        self.assertEqual(len(table), 8)
        self.assertEqual(table["patient_id"].nunique(), 7)
        # The treatment column of the layout is not carried over.
        self.assertEqual(tuple(table.columns), SEER_COLUMNS)

    def test_year_filter(self):
        table = self.prepare(year_range=COHORT)
        self.assertNotIn("00000006", set(table["patient_id"]))
        self.assertIn("00000007", set(table["patient_id"]))

    def test_site_filter(self):
        table = self.prepare(site_recodes=LUNG)
        self.assertIn("00000006", set(table["patient_id"]))
        self.assertNotIn("00000007", set(table["patient_id"]))

    def test_codes_and_blanks(self):
        table = self.prepare(cases_path=CASES).set_index("case_id")
        # Leading zeros of the patient ID are kept; codes become integers.
        self.assertEqual(table.loc["00000002_1", "patient_id"], "00000002")
        self.assertEqual(table.loc["00000002_1", "survival_months"], "12")
        self.assertEqual(table.loc["00000002_1", "registry"], "1502")
        # Blanks stay empty; special codes are kept as they are.
        self.assertEqual(table.loc["00000003_1", "state_county"], "")
        self.assertEqual(table.loc["00000004_1", "grade"], "")
        self.assertEqual(table.loc["00000005_1", "survival_months"], "9999")
        self.assertEqual(table.loc["00000005_1", "age_at_diagnosis"], "999")

    def test_incidence_without_layout_raises(self):
        with self.assertRaises(ValueError):
            prepare_seer(self.tmp, incidence_path=INCIDENCE)

    def test_missing_incidence_raises(self):
        with self.assertRaises(FileNotFoundError):
            self.prepare(incidence_path=self.tmp / "missing")
        with self.assertRaises(FileNotFoundError):
            self.prepare(incidence_path=self.tmp)  # no .TXT files

    def test_missing_required_variable_raises(self):
        layout = self.tmp / "layout.sas"
        layout.write_text(
            "\n".join(
                line
                for line in SAS.read_text().splitlines()
                if "SRV_TIME_MON " not in line
            )
        )
        with self.assertRaisesRegex(ValueError, "survival_months"):
            self.prepare(sas_path=layout)


class TestPrepareFromExport(TempDirTestCase):
    def test_export_without_optional_column(self):
        export = labelled_export(self.tmp, drop=("state_county",))
        with self.assertLogs("pyhealth.datasets.seer", level="WARNING"):
            path = prepare_seer(self.tmp / "out", export_path=export)
        table = read_str(path)
        expected = read_str(FIXTURE / "seer.csv")
        expected["state_county"] = ""
        pd.testing.assert_frame_equal(table, expected)

    def test_tab_separated_export(self):
        table = read_str(labelled_export(self.tmp))
        path = self.tmp / "export.txt"
        table.to_csv(path, sep="\t", index=False)
        out = read_str(prepare_seer(self.tmp / "out", export_path=path))
        pd.testing.assert_frame_equal(out, read_str(FIXTURE / "seer.csv"))

    def test_labels_instead_of_codes_raise(self):
        table = read_str(labelled_export(self.tmp))
        table["SEER cause-specific death classification"] = (
            "Alive or dead of other cause"
        )
        table.to_csv(self.tmp / "export.csv", index=False)
        with self.assertRaisesRegex(ValueError, "codes"):
            prepare_seer(self.tmp / "out", export_path=self.tmp / "export.csv")

    def test_duplicate_case_raises(self):
        table = read_str(labelled_export(self.tmp))
        pd.concat([table, table.iloc[[0]]]).to_csv(self.tmp / "export.csv", index=False)
        with self.assertRaisesRegex(ValueError, "once"):
            prepare_seer(self.tmp / "out", export_path=self.tmp / "export.csv")

    def test_negative_survival_raises(self):
        table = read_str(labelled_export(self.tmp))
        table.loc[0, "Survival months"] = "-1"
        table.to_csv(self.tmp / "export.csv", index=False)
        with self.assertRaisesRegex(ValueError, "negative"):
            prepare_seer(self.tmp / "out", export_path=self.tmp / "export.csv")

    def test_blank_patient_id_raises(self):
        table = read_str(labelled_export(self.tmp))
        table.loc[0, "Patient ID"] = ""
        table.to_csv(self.tmp / "export.csv", index=False)
        with self.assertRaisesRegex(ValueError, "patient_id"):
            prepare_seer(self.tmp / "out", export_path=self.tmp / "export.csv")

    def test_case_list_needs_both_keys(self):
        cases = self.tmp / "cases.csv"
        pd.DataFrame({"Patient ID": ["00000001"]}).to_csv(cases, index=False)
        with self.assertRaises(ValueError):
            prepare_seer(
                self.tmp / "out",
                export_path=labelled_export(self.tmp),
                cases_path=cases,
            )

    def test_source_arguments(self):
        with self.assertRaises(ValueError):
            prepare_seer(self.tmp)
        with self.assertRaises(ValueError):
            prepare_seer(
                self.tmp,
                export_path=labelled_export(self.tmp),
                incidence_path=INCIDENCE,
                sas_path=SAS,
            )


class TestSEERDataset(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(tempfile.mkdtemp())
        shutil.copy(FIXTURE / "seer.csv", cls.root / "seer.csv")
        cls.cache = Path(tempfile.mkdtemp())
        cls.dataset = SEERDataset(root=str(cls.root), cache_dir=str(cls.cache))
        # Load the event table once so each test runs quickly.
        cls.dataset.get_patient("00000001")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root)
        shutil.rmtree(cls.cache)

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_patients(self):
        self.assertEqual(self.dataset.dataset_name, "seer")
        self.assertEqual(
            sorted(self.dataset.unique_patient_ids),
            [f"0000000{i}" for i in range(1, 6)],
        )

    def test_one_event_per_tumor(self):
        events = self.dataset.get_patient("00000001").get_events(event_type="seer")
        self.assertEqual(
            sorted(e["case_id"] for e in events), ["00000001_1", "00000001_2"]
        )
        by_case = {e["case_id"]: e for e in events}
        self.assertEqual(int(by_case["00000001_2"]["cause_of_death_class"]), 1)
        self.assertEqual(int(by_case["00000001_2"]["survival_months"]), 5)

    def test_event_values(self):
        event = self.dataset.get_patient("00000002").get_events(event_type="seer")[0]
        self.assertEqual(int(event["survival_months"]), 12)
        self.assertEqual(int(event["year_of_diagnosis"]), 2006)
        self.assertEqual(int(event["registry"]), 1502)
        blank = self.dataset.get_patient("00000003").get_events(event_type="seer")[0]
        self.assertIn(blank["state_county"], (None, ""))

    def test_missing_csv_without_source_raises(self):
        with self.assertRaisesRegex(FileNotFoundError, "seer.cancer.gov"):
            SEERDataset(root=str(self.tmp))

    def test_missing_column_raises(self):
        read_str(FIXTURE / "seer.csv").drop(columns=["grade"]).to_csv(
            self.tmp / "seer.csv", index=False
        )
        with self.assertRaisesRegex(ValueError, "grade"):
            SEERDataset(root=str(self.tmp))

    def test_builds_csv_from_export(self):
        export = labelled_export(self.tmp)
        root = self.tmp / "root"
        cache = self.tmp / "cache"
        dataset = SEERDataset(
            root=str(root), export_path=str(export), cache_dir=str(cache)
        )
        self.assertTrue((root / "seer.csv").exists())
        self.assertEqual(len(dataset.unique_patient_ids), 5)


if __name__ == "__main__":
    unittest.main()
