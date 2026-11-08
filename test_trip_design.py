"""Integration and rejection tests for fixed-position trip-design workflows.

Run with ``python -m unittest test_trip_design -v``.  All generated fixtures
live in a temporary directory beneath this file's directory and are removed
after the suite.  Synthetic Cp/profile agreement checks consistency only, not
independent physical validation of the shared boundary-layer approximation.
"""

from __future__ import annotations

import contextlib
import copy
import csv
import io
import json
import math
from pathlib import Path
import tempfile
import unittest

from generate_demo import generate
from aero import profile_metrics, roughness_re
from trip_design import main, run_design, select_common_height, validate_config, write_outputs


class TripDesignIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = Path(__file__).resolve().parent
        cls.temporary = tempfile.TemporaryDirectory(prefix="trip-tests-", dir=cls.workspace)
        cls.directory = Path(cls.temporary.name).resolve()
        # Ensure recursive TemporaryDirectory cleanup stays in the test workspace.
        if cls.directory.parent != cls.workspace:
            raise RuntimeError("Temporary test directory escaped the workspace")
        generate(cls.directory)
        cls.config = json.loads((cls.directory / "config.json").read_text(encoding="utf-8"))
        cls.baseline = run_design(cls.directory / "config.json")
        cls.profiles = run_design(cls.directory / "config_profiles.json")

    @classmethod
    def tearDownClass(cls):
        if cls.directory.resolve().parent != cls.workspace:
            raise RuntimeError("Refusing cleanup outside the test workspace")
        cls.temporary.cleanup()

    def config_file(self, cfg, name="modified.json"):
        path = self.directory / (self._testMethodName + "_" + name)
        path.write_text(json.dumps(cfg), encoding="utf-8")
        return path

    def edited_csv(self, source, edit):
        original = self.directory / source
        with original.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            rows = list(reader)
        edit(rows)
        fields = list(dict.fromkeys([*fields, *(column for row in rows for column in row)]))
        output = self.directory / (self._testMethodName + "_" + source)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        return output.name

    def test_synthetic_cp_and_profile_routes_agree(self):
        _, cp_summary, cp_rows = self.baseline
        _, profile_summary, profile_rows = self.profiles
        self.assertEqual(cp_summary["status"], "feasible")
        self.assertEqual(profile_summary["status"], "feasible")
        self.assertEqual(len(cp_rows), 60)
        self.assertEqual(len(cp_summary["design"]), 10)
        for cp, profile in zip(cp_rows, profile_rows):
            self.assertEqual((cp["condition_id"], cp["eta"], cp["surface"]),
                             (profile["condition_id"], profile["eta"], profile["surface"]))
            self.assertLess(abs(cp["k_critical_mm"] / profile["k_critical_mm"] - 1), 1e-5)
            self.assertEqual(cp["selected_height_mm"], profile["selected_height_mm"])
            self.assertGreaterEqual(cp["selected_re_k"], 600)
            self.assertLessEqual(cp["selected_k_delta99"], .8)

    def test_fixed_ten_percent_position_and_swept_coordinates(self):
        cfg, summary, _ = self.baseline
        for row in summary["design"]:
            self.assertEqual(row["x_c"], .1)
            y = row["eta"] * cfg["wing"]["semispan_m"]
            expected_x = y * math.tan(math.radians(30)) + .1 * row["chord_m"]
            self.assertAlmostEqual(row["y_span_m"], y, places=12)
            self.assertAlmostEqual(row["x_global_m"], expected_x, places=12)

    def test_transonic_cp_fallback_rejected_without_heights(self):
        _, summary, rows = run_design(self.directory / "config_transonic_rejected.json")
        self.assertEqual(summary["status"], "not_feasible")
        self.assertTrue(all(r["status"] == "invalid" for r in rows))
        self.assertTrue(all("Mach 0.30" in r["reason"] for r in rows))
        self.assertTrue(all(r["height_mm"] is None for r in summary["design"]))
        self.assertEqual(summary["uniform_height"]["status"], "unavailable")

    def test_missing_condition_never_dropped_from_envelope(self):
        cfg = copy.deepcopy(self.config)
        condition = dict(cfg["conditions"][0], id="missing_condition", alpha_deg=9.)
        cfg["conditions"].append(condition)
        _, summary, rows = run_design(self.config_file(cfg))
        missing = [r for r in rows if r["condition_id"] == "missing_condition"]
        self.assertEqual(len(missing), 10)
        self.assertTrue(all(r["status"] == "invalid" and "missing data" in r["reason"] for r in missing))
        self.assertEqual(summary["status"], "not_feasible")
        self.assertTrue(all(r["height_mm"] is None for r in summary["design"]))

    def test_thin_height_cap_returns_infeasible(self):
        cfg = copy.deepcopy(self.config)
        cfg["trip"]["max_k_delta99"] = .1
        _, summary, rows = run_design(self.config_file(cfg))
        self.assertEqual(summary["status"], "not_feasible")
        self.assertTrue(all(r["status"] == "infeasible" for r in rows))
        self.assertTrue(all("cap" in r["reason"] for r in rows))
        self.assertTrue(all(r["height_mm"] is None for r in summary["design"]))

    def test_common_height_rechecked_against_every_condition(self):
        cfg = copy.deepcopy(self.config)
        cfg["trip"]["max_k_delta99"] = .5
        cfg["trip"]["available_heights_mm"] = []
        _, summary, rows = run_design(self.config_file(cfg))
        # Each case can be tripped separately, but their common height violates
        # the cap in a thinner-layer case: taking only max(k_required) is unsafe.
        self.assertTrue(all(r["status"] == "feasible" for r in rows))
        self.assertEqual(summary["status"], "not_feasible")
        self.assertTrue(any("common" in r["reason"] and "cap" in r["reason"] for r in summary["design"]))

    def test_hardware_rounding_rechecked_against_envelope(self):
        cfg = copy.deepcopy(self.config)
        cfg["trip"]["max_k_delta99"] = .7
        cfg["trip"]["available_heights_mm"] = []
        _, continuous, _ = run_design(self.config_file(cfg, "continuous.json"))
        self.assertEqual(continuous["status"], "feasible")
        cfg["trip"]["available_heights_mm"] = [.25]
        _, rounded, rows = run_design(self.config_file(cfg, "rounded.json"))
        self.assertTrue(all(r["status"] == "feasible" for r in rows))
        self.assertEqual(rounded["status"], "not_feasible")
        rejected = [r for r in rounded["design"] if r["eta"] == .9]
        self.assertTrue(all(r["height_mm"] is None and "common" in r["reason"] for r in rejected))

    def test_nonmonotone_reynolds_searches_later_hardware_and_continuous_roots(self):
        y = [0, .0005, .001, .0015, .002, .003]
        u = [0, 5, 10, 15, 20, 20]
        # This deliberately extreme viscosity variation creates two disjoint
        # admissible Re_k regions and tests the numerical search, not air physics.
        variable = profile_metrics(y, u, [1e-5, 1e-6, 1e-4, 1e-4, 1e-5, 1e-5], 20)
        constant = profile_metrics(y, u, 1e-5, 20)
        self.assertGreater(roughness_re(variable, .0005), 600)
        self.assertLess(roughness_re(variable, .001), 600)
        hardware = select_common_height([variable, constant], .001, [1., 1.5, 1.8], 1., 600)
        self.assertAlmostEqual(hardware, .0018, places=12)
        continuous = select_common_height([variable, constant], .001, [], 1., 600)
        # Between 1.5 and 2 mm: u=10000*k, nu=.00037-.18*k.
        analytical_crossing = (-108 + math.sqrt(108**2 + 4 * 10000 * .222)) / 20000
        self.assertAlmostEqual(continuous, analytical_crossing, places=12)
        self.assertAlmostEqual(roughness_re(variable, continuous), 600, places=7)

    def test_high_freestream_mach_rejected_even_with_slow_local_edge(self):
        cfg = copy.deepcopy(self.config)
        for condition in cfg["conditions"]:
            condition["mach"] = .5
        cp_slow_edge = 2 / (1.4 * .5**2) * (((1 + .2 * .5**2) / (1 + .2 * .15**2))**3.5 - 1)
        def slow_edge(rows):
            for row in rows:
                row["cp"] = str(cp_slow_edge)
        cfg["analysis"]["cp_file"] = self.edited_csv("cp.csv", slow_edge)
        _, _, rows = run_design(self.config_file(cfg))
        self.assertTrue(all(r["status"] == "invalid" and "freestream Mach" in r["reason"] for r in rows))

    def test_non_grit_types_require_explicit_calibration(self):
        for trip_type in ("zigzag", "wire", "custom"):
            cfg = copy.deepcopy(self.config)
            cfg["trip"]["type"] = trip_type
            with self.subTest(trip_type=trip_type), self.assertRaisesRegex(ValueError, "require trip.calibration"):
                validate_config(cfg)
            cfg["trip"]["calibration"] = {
                "critical_re_k": 700., "min_re_x": 1e5, "max_edge_mach": .5,
                "notes": "Test-only explicit calibration; not an empirical recommendation",
            }
            self.assertEqual(validate_config(cfg)["critical_re_k"], 700.)

    def test_downstream_cp_start_requires_initial_thickness(self):
        cfg = copy.deepcopy(self.config)
        def clear_initial(rows):
            for row in rows:
                row["theta0_m"] = ""
        cfg["analysis"]["cp_file"] = self.edited_csv("cp.csv", clear_initial)
        _, summary, rows = run_design(self.config_file(cfg))
        self.assertEqual(summary["status"], "not_feasible")
        self.assertTrue(all(r["status"] == "invalid" and "theta0_m" in r["reason"] for r in rows))

    def test_cp_and_profiles_do_not_extrapolate_position(self):
        cfg = copy.deepcopy(self.config)
        cfg["trip"]["position_xc"] = {"upper": .4, "lower": .4}
        _, _, rows = run_design(self.config_file(cfg, "cp.json"))
        self.assertTrue(all(r["status"] == "invalid" and "no extrapolation" in r["reason"] for r in rows))
        cfg["analysis"] = {"source": "profiles", "profiles_file": "profiles.csv"}
        cfg["trip"]["position_xc"] = {"upper": .11, "lower": .11}
        _, _, rows = run_design(self.config_file(cfg, "profiles.json"))
        self.assertTrue(all(r["status"] == "invalid" and "no supplied profile" in r["reason"] for r in rows))

    def test_false_laminar_confirmation_invalidates_local_design(self):
        cfg = copy.deepcopy(self.config)
        def unconfirm_first(rows):
            rows[0]["laminar_confirmed"] = "false"
        cfg["analysis"] = {
            "source": "profiles",
            "profiles_file": self.edited_csv("profiles.csv", unconfirm_first),
        }
        _, summary, rows = run_design(self.config_file(cfg))
        invalid = [r for r in rows if r["status"] == "invalid"]
        self.assertEqual(len(invalid), 1)
        self.assertIn("confirm undisturbed laminar", invalid[0]["reason"])
        self.assertEqual(summary["status"], "not_feasible")
        local = [r for r in summary["design"] if r["eta"] == .2 and r["surface"] == "upper"]
        self.assertIsNone(local[0]["height_mm"])

    def test_known_transition_extent_rejected(self):
        cfg = copy.deepcopy(self.config)
        cfg["analysis"]["laminar_until_xc"] = .09
        _, _, rows = run_design(self.config_file(cfg))
        self.assertTrue(all(r["status"] == "invalid" and "laminar extent" in r["reason"] for r in rows))

    def test_nonfinite_laminar_limits_rejected(self):
        for invalid in (float("nan"), float("inf"), -0.1, 1.1):
            cfg = copy.deepcopy(self.config)
            cfg["analysis"]["laminar_until_xc"] = invalid
            with self.subTest(limit=invalid), self.assertRaisesRegex(ValueError, "laminar_until_xc"):
                validate_config(cfg)
        cfg = copy.deepcopy(self.config)
        def corrupt_limit(rows):
            # A non-first row must be checked too; otherwise malformed upstream
            # transition metadata could silently disappear.
            rows[5]["known_transition_xc"] = "nan"
        cfg["analysis"]["cp_file"] = self.edited_csv("cp.csv", corrupt_limit)
        _, summary, rows = run_design(self.config_file(cfg))
        invalid_rows = [r for r in rows if r["status"] == "invalid"]
        self.assertEqual(len(invalid_rows), 1)
        self.assertIn("known_transition_xc", invalid_rows[0]["reason"])
        self.assertEqual(summary["status"], "not_feasible")

    def test_inconsistent_known_transition_metadata_rejected(self):
        cfg = copy.deepcopy(self.config)
        def conflicting_limits(rows):
            rows[0]["known_transition_xc"] = ".2"
            rows[5]["known_transition_xc"] = ".3"
        cfg["analysis"]["cp_file"] = self.edited_csv("cp.csv", conflicting_limits)
        _, _, rows = run_design(self.config_file(cfg))
        invalid_rows = [r for r in rows if r["status"] == "invalid"]
        self.assertEqual(len(invalid_rows), 1)
        self.assertIn("must be consistent", invalid_rows[0]["reason"])

    def test_low_running_reynolds_not_used_with_default_grit_rule(self):
        cfg = copy.deepcopy(self.config)
        for condition in cfg["conditions"]:
            condition["reynolds"] = 1e4
        _, _, rows = run_design(self.config_file(cfg))
        self.assertTrue(all(r["status"] == "invalid" and "below calibration minimum" in r["reason"] for r in rows))

    def test_output_serialization_and_negative_case_exit_status(self):
        cfg, summary, rows = copy.deepcopy(self.baseline)
        output = self.directory / "written_results"
        write_outputs(output, cfg, summary, rows, plot=False)
        for filename in ("design.csv", "conditions.csv", "summary.json", "report.html"):
            self.assertTrue((output / filename).is_file())
        serialized = json.loads((output / "summary.json").read_text(encoding="utf-8"))
        self.assertEqual(serialized["status"], "feasible")
        self.assertFalse(serialized["plot_created"])
        self.assertEqual(len(serialized["inputs"][0]["sha256"]), 64)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = main(["--config", str(self.directory / "config_transonic_rejected.json"),
                         "--output", str(self.directory / "negative_results"), "--no-plot"])
        self.assertEqual(code, 2)
        self.assertTrue((self.directory / "negative_results" / "summary.json").is_file())


if __name__ == "__main__":
    unittest.main()
