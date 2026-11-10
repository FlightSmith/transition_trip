"""Report-only tests from synthetic saved BL integrals; no XFOIL execution."""

from __future__ import annotations

import copy
import csv
import html
import json
from pathlib import Path
import shutil
import sys
import tempfile
import types
import unittest
from unittest import mock

from aero import tunnel_state
from scaled_wing_report import (DEFAULT_AIRFOIL_NAME, POSITIONS, REYNOLDS,
                                _plot_comparison, _validate_inputs, summarize)


class ScaledWingReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = Path(__file__).resolve().parent
        cls.temporary = tempfile.TemporaryDirectory(prefix="scaled-report-tests-", dir=cls.workspace)
        cls.directory = Path(cls.temporary.name).resolve()
        if cls.directory.parent != cls.workspace:
            raise RuntimeError("Temporary test directory escaped workspace")
        source = cls.directory / "all_conditions"
        source.mkdir()
        conditions = [{"id": f"re{int(re/1e6)}_a{alpha}", "reynolds": re,
                       "alpha_deg": alpha, "mach": .2, "temperature_k": 288.15}
                      for re in REYNOLDS for alpha in (0, 1, 2)]
        stations = [{"eta": i / 20, "chord_m": .695 - .57 * i / 20,
                     "y_m": -.166625 - .614625 * i / 20,
                     "x_le_m": 1.16 + .409625 * i / 20}
                    for i in range(21)]
        cls.configs = {}
        for xc in POSITIONS:
            cfg = {"schema_version": 1, "name": "Synthetic report-test fixture",
                   "data_origin": "Test fixture only: constant integral quantities, not aerodynamic data",
                   "wing": {"semispan_m": .614625, "reference_length_m": .578625,
                            "sweep_le_deg": 33.68199959181396, "stations": stations},
                   "conditions": conditions, "trip": {"type": "grit", "height_margin": 1.1,
                   "max_k_delta99": .8, "available_heights_mm": [],
                   "position_xc": {"upper": xc, "lower": xc}},
                   "analysis": {"source": "integral", "integral_file": "boundary_layer.csv"}}
            cls.configs[xc] = cfg
            (source / f"config_xc_{xc:g}.json").write_text(json.dumps(cfg), encoding="utf-8")
        cls.records = []
        for condition in conditions:
            flow = tunnel_state(.2, condition["reynolds"], .578625)
            for station in stations:
                for xc in POSITIONS:
                    for surface in ("upper", "lower"):
                        cls.records.append({"condition_id": condition["id"], "eta": station["eta"],
                            "surface": surface, "x_c": xc, "status": "valid", "reason": "",
                            "theta_m": .00006, "shape_factor": 2.6, "nu_m2_s": flow["nu_inf_m2_s"],
                            "ue_m_s": flow["U_inf_m_s"], "edge_mach": .2,
                            "laminar_confirmed": "true", "transition_xc": .8, "cf": .001})
        with (source / "boundary_layer.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(cls.records[0]))
            writer.writeheader()
            writer.writerows(cls.records)
        (source / "workflow_summary.json").write_text(json.dumps({"diagnostics": [], "solver_plan": []}), encoding="utf-8")
        with mock.patch("xfoil_backend.run_xfoil", side_effect=AssertionError("report must not execute XFOIL")) as runner:
            cls.result = summarize(cls.directory, plot=False)
            runner.assert_not_called()
        with (cls.directory / "comparison.csv").open(newline="", encoding="utf-8") as handle:
            cls.comparison = list(csv.DictReader(handle))

    @classmethod
    def tearDownClass(cls):
        if cls.directory.resolve().parent != cls.workspace:
            raise RuntimeError("Refusing cleanup outside test workspace")
        cls.temporary.cleanup()

    def test_per_re_data_and_conditions_are_filtered_together(self):
        self.assertEqual(self.result["comparison_rows"], 24)
        for re in REYNOLDS:
            scope = f"re{int(re/1e6)}e6"
            directory = self.directory / scope
            cfg = json.loads((directory / "config_xc_0.1.json").read_text())
            self.assertEqual(len(cfg["conditions"]), 3)
            self.assertTrue(all(c["reynolds"] == re for c in cfg["conditions"]))
            ids = {c["id"] for c in cfg["conditions"]}
            with (directory / "boundary_layer.csv").open(newline="", encoding="utf-8") as handle:
                records = list(csv.DictReader(handle))
            self.assertEqual(len(records), 3 * 21 * 3 * 2)
            self.assertEqual({r["condition_id"] for r in records}, ids)
            self.assertTrue((directory / "xc_0.1" / "report.html").is_file())

    def test_calibration_masks_produce_expected_partial_span_counts(self):
        counts = {"re3e6": {.07: 15, .1: 18, .12: 19},
                  "re5e6": {.07: 19, .1: 21, .12: 21},
                  "re8e6": {.07: 21, .1: 21, .12: 21},
                  "all_conditions": {.07: 15, .1: 18, .12: 19}}
        for row in self.comparison:
            self.assertEqual(int(row["valid_station_count"]), counts[row["scope"]][float(row["position_xc"])])
            self.assertEqual(int(row["total_station_count"]), 21)

    def test_partial_sets_never_get_uniform_height_or_numeric_invalid_tip(self):
        for row in self.comparison:
            if int(row["valid_station_count"]) < 21:
                self.assertEqual(row["uniform_all_stations_height_mm"], "")
                self.assertEqual(row["tip_height_mm"], "")
                self.assertEqual(row["tip_status"], "invalid")
                self.assertIn("below calibration minimum", row["tip_reason"])
            else:
                self.assertGreater(float(row["uniform_all_stations_height_mm"]), 0.)
                self.assertEqual(row["uniform_status"], "feasible")
            self.assertGreater(float(row["min_valid_height_mm"]), 0.)
            self.assertGreaterEqual(float(row["max_valid_height_mm"]), float(row["min_valid_height_mm"]))

    def test_report_and_audit_explain_partial_coverage_and_saved_data(self):
        html = (self.directory / "report.html").read_text(encoding="utf-8")
        self.assertIn("Partial coverage is not a full-span design", html)
        self.assertIn("does not rerun XFOIL", html)
        audit = json.loads((self.directory / "comparison_summary.json").read_text())
        self.assertFalse(audit["chart_created"])
        self.assertEqual(len(audit["groups"]), 12)
        self.assertEqual(len(audit["inputs"]), 5)
        self.assertTrue(all(not Path(item["path"]).is_absolute() for item in audit["inputs"]))

    def test_legacy_inputs_default_to_naca_name(self):
        self.assertEqual(self.result["airfoil_name"], "NACA 64A204")
        report = (self.directory / "report.html").read_text(encoding="utf-8")
        self.assertIn("<h1>1:8 NACA 64A204 swept-wing trip comparison</h1>", report)
        audit = json.loads((self.directory / "comparison_summary.json").read_text())
        self.assertEqual(audit["airfoil_name"], DEFAULT_AIRFOIL_NAME)
        for row in self.comparison:
            self.assertEqual(row["airfoil_name"], DEFAULT_AIRFOIL_NAME)
        cfg = json.loads((self.directory / "re3e6" / "config_xc_0.1.json").read_text())
        self.assertEqual(cfg["airfoil_name"], DEFAULT_AIRFOIL_NAME)
        self.assertIn(DEFAULT_AIRFOIL_NAME, cfg["name"])

    def test_custom_airfoil_name_propagates_and_html_escapes(self):
        airfoil_name = 'Custom <b>foil</b> & "quoted" $x$'
        output = self.directory / "custom_airfoil"
        source = output / "all_conditions"
        source.mkdir(parents=True)
        for position, original in self.configs.items():
            cfg = copy.deepcopy(original)
            cfg["airfoil_name"] = airfoil_name
            cfg["name"] = "Historical NACA 64A204 display name"
            (source / f"config_xc_{position:g}.json").write_text(json.dumps(cfg), encoding="utf-8")
        for name in ("boundary_layer.csv", "workflow_summary.json"):
            shutil.copyfile(self.directory / "all_conditions" / name, source / name)
        with mock.patch("xfoil_backend.run_xfoil", side_effect=AssertionError("report must not execute XFOIL")):
            result = summarize(output, plot=False)
        self.assertEqual(result["airfoil_name"], airfoil_name)
        escaped = html.escape(airfoil_name)
        report = (output / "report.html").read_text(encoding="utf-8")
        self.assertIn(f"<title>1:8 {escaped} swept-wing trip comparison</title>", report)
        self.assertIn(f"<h1>1:8 {escaped} swept-wing trip comparison</h1>", report)
        self.assertNotIn(airfoil_name, report)
        for scope in ("all_conditions", "re3e6", "re5e6", "re8e6"):
            summary = json.loads((output / scope / "xc_0.1" / "summary.json").read_text())
            self.assertIn(airfoil_name, summary["name"])
            self.assertNotIn("NACA", summary["name"])
            individual = (output / scope / "xc_0.1" / "report.html").read_text(encoding="utf-8")
            self.assertIn(escaped, individual)
            self.assertNotIn(airfoil_name, individual)
            if scope != "all_conditions":
                derived = json.loads((output / scope / "config_xc_0.1.json").read_text())
                self.assertEqual(derived["airfoil_name"], airfoil_name)
                self.assertIn(airfoil_name, derived["name"])
        audit = json.loads((output / "comparison_summary.json").read_text())
        self.assertEqual(audit["airfoil_name"], airfoil_name)
        with (output / "comparison.csv").open(newline="", encoding="utf-8") as handle:
            self.assertTrue(all(row["airfoil_name"] == airfoil_name for row in csv.DictReader(handle)))

    def test_plot_uses_literal_custom_airfoil_name(self):
        summaries = {}
        for scope in ("re3e6", "re5e6", "re8e6"):
            for xc in POSITIONS:
                summaries[(scope, xc)] = json.loads((self.directory / scope / f"xc_{xc:g}" / "summary.json").read_text())
        figure = mock.Mock()
        axes = [mock.Mock() for _ in POSITIONS]
        for axis in axes:
            axis.get_legend_handles_labels.return_value = ([], [])
        plotting = types.ModuleType("matplotlib.pyplot")
        plotting.rcParams = {}
        plotting.subplots = mock.Mock(return_value=(figure, axes))
        plotting.close = mock.Mock()
        matplotlib = types.ModuleType("matplotlib")
        matplotlib.use = mock.Mock()
        matplotlib.pyplot = plotting
        name = "Custom $literal$ airfoil"
        with mock.patch.dict(sys.modules, {"matplotlib": matplotlib, "matplotlib.pyplot": plotting}), \
                mock.patch.dict("os.environ", {}):
            self.assertTrue(_plot_comparison(self.directory, summaries, name))
        self.assertEqual(figure.suptitle.call_args.args[0], f"{name} wing | 1:8 model | Mach 0.20")
        self.assertFalse(figure.suptitle.call_args.kwargs["parse_math"])
        self.assertEqual(figure.savefig.call_count, 2)

    def test_airfoil_names_must_be_nonempty_and_consistent(self):
        configs = copy.deepcopy(self.configs)
        configs[.1]["airfoil_name"] = "Different airfoil"
        with self.assertRaisesRegex(ValueError, "same airfoil_name"):
            _validate_inputs(configs, self.records)
        for invalid in ("", "  ", None, 123, ["airfoil"]):
            configs = copy.deepcopy(self.configs)
            configs[.07]["airfoil_name"] = invalid
            with self.subTest(name=invalid), self.assertRaisesRegex(ValueError, "nonempty string"):
                _validate_inputs(configs, self.records)

    def test_missing_condition_matrix_and_unknown_saved_ids_rejected(self):
        configs = copy.deepcopy(self.configs)
        configs[.07]["conditions"].pop()
        with self.assertRaisesRegex(ValueError, "nine unique"):
            _validate_inputs(configs, self.records)
        records = [*self.records, dict(self.records[0], condition_id="unknown")]
        with self.assertRaisesRegex(ValueError, "unknown condition"):
            _validate_inputs(self.configs, records)


if __name__ == "__main__":
    unittest.main()
