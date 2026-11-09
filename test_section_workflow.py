"""Integration checks for raw section/Cp import and optional solver orchestration.

XFOIL is mocked in every solver test; this suite never launches an executable.
Raw fixtures are copies of the user's supplied synthetic_Cp files.
"""

from __future__ import annotations

import contextlib
import copy
import csv
import io
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

from aero import edge_state, tunnel_state
from curve_io import load_section, map_curve_to_section, read_curve
from section_workflow import (default_settings, discover, execute,
                              initialized_cp_rows, validate_settings, write_json)


class SectionWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = Path(__file__).resolve().parent
        cls.raw_directory = cls.workspace / "synthetic_Cp"
        if not cls.raw_directory.is_dir():
            raise unittest.SkipTest("synthetic_Cp raw example directory is unavailable")

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="section-workflow-tests-", dir=self.workspace)
        self.directory = Path(self.temporary.name).resolve()
        self.assertEqual(self.directory.parent, self.workspace)
        self.inputs = self.copy_inputs("inputs")
        self.settings = default_settings("inputs")

    def tearDown(self):
        if self.directory.resolve().parent != self.workspace:
            raise RuntimeError("Refusing cleanup outside test workspace")
        self.temporary.cleanup()

    def copy_inputs(self, name):
        destination = self.directory / name
        destination.mkdir()
        for source in self.raw_directory.iterdir():
            if source.is_file():
                shutil.copyfile(source, destination / source.name)
        return destination

    def run_workflow(self, settings=None, output_name="output", **options):
        path = self.directory / (output_name + "_settings.json")
        write_json(path, settings or self.settings)
        output = self.directory / output_name
        with contextlib.redirect_stdout(io.StringIO()):
            code = execute(path, output, no_plot=True, **options)
        return code, output

    @staticmethod
    def read_table(path):
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    @staticmethod
    def replace_curve_values(path, points):
        header = path.read_text(encoding="utf-8").splitlines()[:7]
        lines = [*header, *(f"{x:.16E} {y:.16E}" for x, y in points), "*END"]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    def test_raw_curves_use_their_actual_geometry_and_all_states(self):
        code, output = self.run_workflow()
        self.assertEqual(code, 0)
        metadata = json.loads((output / "workflow_summary.json").read_text())
        self.assertEqual(len(metadata["imported_files"]), 6)
        self.assertEqual(len(metadata["diagnostics"]), 12)
        geometry = metadata["geometry"]
        self.assertEqual([g["y_m"] for g in geometry], [-1.25, -5.])
        self.assertEqual([g["chord_m"] for g in geometry], [1.75, 1.])
        cp_rows = self.read_table(output / "cp.csv")
        self.assertEqual(len(cp_rows), 3 * 2 * 2 * 41)
        self.assertEqual({r["condition_id"] for r in cp_rows}, {"state1", "state2", "state3"})
        cfg = json.loads((output / "config_xc_0.1.json").read_text())
        self.assertEqual(cfg["wing"]["semispan_m"], 3.75)
        self.assertEqual(cfg["wing"]["reference_length_m"], 4.629)
        self.assertEqual([c["alpha_deg"] for c in cfg["conditions"]], [0., 1., 2.])
        for label in ("xc_0.07", "xc_0.1", "xc_0.12"):
            result = json.loads((output / label / "summary.json").read_text())
            self.assertEqual(result["status"], "feasible")
            self.assertEqual(len(result["design"]), 4)

    def test_absolute_pressure_import_is_equivalent_to_cp(self):
        _, cp_output = self.run_workflow(output_name="cp")
        pressure_inputs = self.copy_inputs("pressures")
        q, reference = 500., 100000.
        for path in pressure_inputs.glob("cp_dist*"):
            points = read_curve(path)
            points[:, 1] = reference + q * points[:, 1]
            self.replace_curve_values(path, points)
        settings = copy.deepcopy(self.settings)
        settings.update(input_directory="pressures", value_kind="pressure",
                        pressure_reference={"dynamic_pressure_pa": q, "static_pressure_pa": reference})
        code, pressure_output = self.run_workflow(settings, "pressure")
        self.assertEqual(code, 0)
        cp_rows = self.read_table(cp_output / "cp.csv")
        pressure_rows = self.read_table(pressure_output / "cp.csv")
        np.testing.assert_allclose([float(r["cp"]) for r in pressure_rows],
                                   [float(r["cp"]) for r in cp_rows], atol=3e-14, rtol=1e-10)
        for label in ("xc_0.07", "xc_0.1", "xc_0.12"):
            cp_result = json.loads((cp_output / label / "summary.json").read_text())
            pressure_result = json.loads((pressure_output / label / "summary.json").read_text())
            np.testing.assert_allclose([r["height_mm"] for r in cp_result["design"]],
                                       [r["height_mm"] for r in pressure_result["design"]], rtol=1e-9)

    def test_pressure_import_requires_explicit_references(self):
        settings = copy.deepcopy(self.settings)
        settings["value_kind"] = "pressure"
        with self.assertRaisesRegex(ValueError, "requires both"):
            self.run_workflow(settings)
        settings["pressure_reference"] = {"dynamic_pressure_pa": 500.}
        with self.assertRaisesRegex(ValueError, "requires both"):
            self.run_workflow(settings, "missing_static")

    def test_missing_configured_state_and_duplicate_files_rejected(self):
        settings = copy.deepcopy(self.settings)
        settings["states"].append({"state": 4, "alpha_deg": 3.})
        with self.assertRaisesRegex(ValueError, "missing curve for state4"):
            self.run_workflow(settings)
        source = self.inputs / "section_state1_station-1.250"
        shutil.copyfile(source, source.with_name(source.name + ".dat"))
        with self.assertRaisesRegex(ValueError, "duplicate state/station"):
            discover(self.inputs)

    def test_changing_geometry_between_states_rejected(self):
        source = self.inputs / "section_state1_station-1.250"
        different = source.with_name("section_state2_station-1.250")
        shutil.copyfile(source, different)
        points = read_curve(different)
        points[:, 1] += .01
        self.replace_curve_values(different, points)
        with self.assertRaisesRegex(ValueError, "State-varying section geometry"):
            self.run_workflow()

    def test_reference_state_and_data_start_line_require_positive_integers(self):
        for field in ("section_state", "data_start_line"):
            for value in (0, -1, 1.5, float("nan")):
                settings = copy.deepcopy(self.settings)
                settings[field] = value
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    validate_settings(settings)

    def test_custom_data_start_line_applies_to_sections_and_cp(self):
        for path in self.inputs.iterdir():
            lines = path.read_text().splitlines()
            lines.insert(7, "42 0")  # Numeric header metadata, not a contour point.
            path.write_text("\n".join(lines) + "\n")
        settings = copy.deepcopy(self.settings)
        settings["data_start_line"] = 9
        code, output = self.run_workflow(settings)
        self.assertEqual(code, 0)
        metadata = json.loads((output / "workflow_summary.json").read_text())
        self.assertEqual([g["chord_m"] for g in metadata["geometry"]], [1.75, 1.])

    def test_initialization_retains_true_arc_distance_and_state_override(self):
        section = load_section(self.inputs / "section_state1_station-1.250")
        surface = map_curve_to_section(self.inputs / "cp_dist_state1_station-1.250", section)["upper"]
        flow = tunnel_state(.2, 8e6, 4.629)
        settings = copy.deepcopy(self.settings)
        settings["cp_initialization"] = {"method": "provided", "theta0_m": 2e-5}
        state = {"cp_initialization": {"method": "provided", "theta0_m": 3e-5}}
        x, cp, s, theta0, method = initialized_cp_rows(surface, settings, state, flow)
        self.assertEqual(theta0, 3e-5)
        self.assertEqual(method, "provided")
        self.assertEqual(x[0], surface["x_c"][1])
        self.assertEqual(s[0], surface["s_m"][1])
        self.assertGreater(s[0], 0.)
        self.assertGreater(s[0], x[0] * section["chord_m"])
        settings["cp_initialization"] = {"method": "stagnation_ramp"}
        xx, pp, ss, theta, _ = initialized_cp_rows(surface, settings, {}, flow)
        ue = float(edge_state([pp[0]], flow, .2)["U_e_m_s"][0])
        self.assertAlmostEqual(theta, math.sqrt(.075 * flow["nu_inf_m2_s"] * ss[0] / ue))
        self.assertEqual(ss[0], s[0])
        incomplete = {k: v[1:] for k, v in surface.items() if isinstance(v, np.ndarray)}
        with self.assertRaisesRegex(ValueError, "actual geometric leading edge"):
            initialized_cp_rows(incomplete, self.settings, {}, flow)

    def test_prepare_only_xfoil_does_not_call_external_solver(self):
        settings = copy.deepcopy(self.settings)
        settings["mode"] = "xfoil"
        with mock.patch("xfoil_backend.run_xfoil", side_effect=AssertionError("must not run")) as runner:
            code, output = self.run_workflow(settings, prepare_only=True)
        self.assertEqual(code, 0)
        runner.assert_not_called()
        metadata = json.loads((output / "workflow_summary.json").read_text())
        self.assertEqual(metadata["status"], "prepared")
        self.assertEqual(len(metadata["solver_plan"]), 6)
        self.assertEqual(len(list((output / "airfoils").glob("*.dat"))), 2)
        self.assertFalse((output / "boundary_layer.csv").exists())
        for plan in metadata["solver_plan"]:
            chord = 1.75 if plan["y_m"] == -1.25 else 1.
            self.assertAlmostEqual(plan["reynolds_local"], 8e6 * chord / 4.629)
            self.assertFalse(Path(plan["airfoil_file"]).is_absolute())

    def test_global_angle_reference_subtracts_section_incidence(self):
        angle = math.radians(2.)
        rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
        for path in self.inputs.glob("section*"):
            points = read_curve(path)
            leading = points[np.argmin(points[:, 0])].copy()
            rotated = (points - leading) @ rotation.T + leading
            self.replace_curve_values(path, rotated)
        settings = copy.deepcopy(self.settings)
        settings.update(mode="xfoil", alpha_reference="global_x", states=[{"state": 1, "alpha_deg": 5.}])
        with mock.patch("xfoil_backend.run_xfoil", side_effect=AssertionError("must not run")):
            _, output = self.run_workflow(settings, prepare_only=True)
        metadata = json.loads((output / "workflow_summary.json").read_text())
        for plan in metadata["solver_plan"]:
            self.assertAlmostEqual(plan["section_incidence_deg"], 2., places=10)
            self.assertAlmostEqual(plan["alpha_local_deg"], 3., places=10)
        for geometry in metadata["geometry"]:
            self.assertAlmostEqual(geometry["streamwise_chord_m"], geometry["chord_m"] * math.cos(angle))

    def test_nonconverged_xfoil_result_becomes_invalid_reports(self):
        settings = copy.deepcopy(self.settings)
        settings["mode"] = "xfoil"
        failed = {"converged": False, "reason": "deliberate mock convergence failure"}
        with mock.patch("xfoil_backend.run_xfoil", return_value=failed) as runner, \
                mock.patch("xfoil_backend.sample_laminar_station") as sample:
            code, output = self.run_workflow(settings)
        self.assertEqual(code, 2)
        self.assertEqual(runner.call_count, 6)
        sample.assert_not_called()
        metadata = json.loads((output / "workflow_summary.json").read_text())
        self.assertEqual(len(metadata["diagnostics"]), 6)
        self.assertTrue(all(d["status"] == "failed" for d in metadata["diagnostics"]))
        rows = self.read_table(output / "boundary_layer.csv")
        self.assertEqual(len(rows), 36)
        self.assertTrue(all(r["status"] == "invalid" and r["laminar_confirmed"] == "false" for r in rows))
        for label in ("xc_0.07", "xc_0.1", "xc_0.12"):
            summary = json.loads((output / label / "summary.json").read_text())
            self.assertEqual(summary["status"], "not_feasible")
            self.assertTrue(all(r["height_mm"] is None for r in summary["design"]))

    def test_mocked_solver_integrals_scale_by_true_local_chord(self):
        settings = copy.deepcopy(self.settings)
        settings.update(mode="xfoil", positions_xc=[.1], states=[{"state": 1, "alpha_deg": 0.}])
        solved = {"converged": True, "boundary_layer": {"h_column_interpretation": "Hk"},
                  "cp": {"x_c": np.array([1., .5, 0., .5, 1.]), "cp": np.zeros(5)}}
        sample = {"theta_c": .0001, "delta_star_c": .00026, "h_kinematic": 2.6,
                  "h_mass": 2.64, "cf": .001, "ue_over_uinf": 1.1,
                  "edge_mach": .22, "transition_xc": .6}
        with mock.patch("xfoil_backend.run_xfoil", return_value=solved) as runner, \
                mock.patch("xfoil_backend.sample_laminar_station", return_value=sample):
            code, output = self.run_workflow(settings)
        self.assertEqual(code, 0)
        self.assertEqual(runner.call_count, 2)
        rows = self.read_table(output / "boundary_layer.csv")
        self.assertEqual(len(rows), 4)
        for row in rows:
            chord = 1.75 if float(row["eta"]) == 0 else 1.
            self.assertAlmostEqual(float(row["theta_m"]), sample["theta_c"] * chord)
            self.assertAlmostEqual(float(row["delta_star_m"]), sample["delta_star_c"] * chord)
            self.assertEqual(float(row["shape_factor"]), 2.6)
            self.assertEqual(float(row["h_mass"]), 2.64)
            self.assertEqual(row["laminar_confirmed"], "true")


if __name__ == "__main__":
    unittest.main()
