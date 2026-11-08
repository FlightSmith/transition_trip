"""Portable Cp-input example checks; generated Cp is placeholder data only."""

from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from generate_left_wing_cp import CP_FIELDS, generate
from trip_design import run_design


class LeftWingCpExampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = Path(__file__).resolve().parent
        cls.temporary = tempfile.TemporaryDirectory(prefix="left-wing-cp-tests-", dir=cls.workspace)
        cls.directory = Path(cls.temporary.name).resolve()
        if cls.directory.parent != cls.workspace:
            raise RuntimeError("Temporary test directory escaped the workspace")
        cls.generated = generate(cls.directory)
        cls.results = {code: run_design(cls.directory / f"config_xc{code}.json")
                       for code in ("07", "10", "12")}
        with (cls.directory / "cp.csv").open(newline="", encoding="utf-8") as handle:
            cls.cp_rows = list(csv.DictReader(handle))

    @classmethod
    def tearDownClass(cls):
        if cls.directory.resolve().parent != cls.workspace:
            raise RuntimeError("Refusing cleanup outside the test workspace")
        cls.temporary.cleanup()

    def modified_config(self, cp_filename):
        cfg = copy.deepcopy(self.results["10"][0])
        cfg["analysis"]["cp_file"] = cp_filename
        path = self.directory / (self._testMethodName + ".json")
        path.write_text(json.dumps(cfg), encoding="utf-8")
        return path

    def test_all_three_positions_run_with_expected_uniform_height_difference(self):
        for code, (_, summary, rows) in self.results.items():
            with self.subTest(position=code):
                self.assertEqual(summary["source"], "cp")
                self.assertEqual(summary["status"], "feasible")
                self.assertEqual(len(summary["design"]), 42)
                self.assertTrue(all(r["status"] == "feasible" for r in rows))
                self.assertTrue(all(r["height_mm"] > 0 for r in summary["design"]))
        # Cp's integral/quartic approximation and numerical Blasius are distinct
        # models; the narrow 7%-chord uniform-height window changes classification.
        self.assertEqual(self.results["07"][1]["uniform_height"]["status"], "infeasible")
        self.assertIsNone(self.results["07"][1]["uniform_height"]["height_mm"])
        self.assertEqual(self.results["10"][1]["uniform_height"]["status"], "feasible")
        self.assertEqual(self.results["12"][1]["uniform_height"]["status"], "feasible")

    def test_placeholder_provenance_and_relative_paths_are_explicit(self):
        for cfg, summary, _ in self.results.values():
            self.assertIn("PLACEHOLDER Cp=0", cfg["data_origin"])
            self.assertIn("NOT ACTUAL WING Cp", summary["data_origin"])
            self.assertEqual(cfg["analysis"], {"source": "cp", "cp_file": "cp.csv"})
            self.assertEqual(cfg["conditions"][0]["id"], "case_01")
            self.assertIn("physical AoA unspecified", cfg["conditions"][0]["alpha_label"])
            self.assertEqual(cfg["wing"]["reference_length_m"], 4.629)
            self.assertEqual(cfg["wing"]["reference_area_m2"], 49.)
            self.assertEqual(cfg["wing"]["stations"][0]["y_m"], -1.333)
            self.assertEqual(cfg["wing"]["stations"][-1]["y_m"], -6.250)

    def test_cp_grid_starts_at_zero_with_only_first_initial_thickness(self):
        self.assertEqual(self.generated["cp_rows"], 7392)
        self.assertTrue(all(float(r["cp"]) == 0 for r in self.cp_rows))
        first_group = [r for r in self.cp_rows if float(r["eta"]) == 0 and r["surface"] == "upper"]
        positions = [float(r["x_c"]) for r in first_group]
        self.assertEqual(len(positions), 176)
        self.assertEqual(positions[0], 0.)
        self.assertEqual(positions[-1], .35)
        for x in (.07, .10, .12):
            self.assertIn(x, positions)
        self.assertEqual(float(first_group[0]["theta0_m"]), 0.)
        self.assertTrue(all(r["theta0_m"] == "" for r in first_group[1:]))
        for r in first_group:
            self.assertAlmostEqual(float(r["s_m"]), float(r["x_c"]) * 5.56)

    def test_missing_cp_file_and_missing_cp_column_are_rejected(self):
        with self.assertRaises(FileNotFoundError):
            run_design(self.modified_config("absent_cp.csv"))
        bad = self.directory / "no_cp_column.csv"
        with bad.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=[k for k in CP_FIELDS if k != "cp"])
            writer.writeheader()
            writer.writerow({k: v for k, v in self.cp_rows[0].items() if k != "cp"})
        with self.assertRaisesRegex(ValueError, "missing columns.*cp"):
            run_design(self.modified_config(bad.name))

    def test_empty_cp_file_and_missing_upstream_thickness_are_rejected(self):
        empty = self.directory / "empty_cp.csv"
        with empty.open("w", newline="", encoding="utf-8") as handle:
            csv.DictWriter(handle, fieldnames=CP_FIELDS).writeheader()
        with self.assertRaisesRegex(ValueError, "empty data"):
            run_design(self.modified_config(empty.name))
        downstream = self.directory / "downstream_cp.csv"
        with downstream.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=CP_FIELDS)
            writer.writeheader()
            writer.writerows(r for r in self.cp_rows if float(r["x_c"]) > 0)
        _, summary, rows = run_design(self.modified_config(downstream.name))
        self.assertEqual(summary["status"], "not_feasible")
        self.assertTrue(all(r["status"] == "invalid" and "theta0_m" in r["reason"] for r in rows))

    def test_generation_is_deterministic_with_lf_files(self):
        names = ["cp.csv", *self.generated["config_files"]]
        before = {name: (self.directory / name).read_bytes() for name in names}
        generate(self.directory)
        for name in names:
            data = (self.directory / name).read_bytes()
            self.assertEqual(data, before[name])
            self.assertNotIn(b"\r\n", data)


if __name__ == "__main__":
    unittest.main()
