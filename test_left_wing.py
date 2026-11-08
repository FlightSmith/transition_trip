"""Independent Blasius checks and user's left-panel baseline integration tests."""

from __future__ import annotations

import json
import math
from pathlib import Path
import tempfile
import unittest

from aero import profile_metrics, required_height, tunnel_state
from generate_left_wing import (
    BLASIUS_WALL_SECOND_DERIVATIVE, MACH, REFERENCE_MAC_M, REYNOLDS_MAC,
    ROOT_CHORD_M, TIP_CHORD_M, TRIP_POSITIONS, blasius_solution, generate,
    geometry_metadata,
)
from trip_design import run_design


class LeftWingBlasiusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = Path(__file__).resolve().parent
        cls.temporary = tempfile.TemporaryDirectory(prefix="left-wing-tests-", dir=cls.workspace)
        cls.directory = Path(cls.temporary.name).resolve()
        if cls.directory.parent != cls.workspace:
            raise RuntimeError("Temporary test directory escaped the workspace")
        cls.generated = generate(cls.directory)
        cls.solution = blasius_solution()
        cls.results = {xc: run_design(cls.directory / f"config_xc{round(xc*100):02d}.json")
                       for xc in TRIP_POSITIONS}

    @classmethod
    def tearDownClass(cls):
        if cls.directory.resolve().parent != cls.workspace:
            raise RuntimeError("Refusing cleanup outside the test workspace")
        cls.temporary.cleanup()

    def test_blasius_integrals_match_classical_values(self):
        solution = self.solution
        # The momentum integral coefficient is exactly 2*f''(0), independently
        # of the numerical quadrature of the sampled velocity profile.
        self.assertLess(abs(solution["theta_coefficient"] - 2 * BLASIUS_WALL_SECOND_DERIVATIVE), 3e-6)
        self.assertAlmostEqual(solution["displacement_coefficient"], 1.72078766, delta=4e-6)
        self.assertAlmostEqual(solution["eta_delta99"], 4.90999, delta=2e-5)
        self.assertLess(abs(solution["raw_edge_velocity_fraction"] - 1), 1e-8)
        self.assertEqual(solution["u_over_U"][0], 0)
        self.assertEqual(solution["u_over_U"][-1], 1)

    def test_profile_sampling_changes_critical_height_by_less_than_ten_ppm(self):
        dense = blasius_solution(sample_step=.001)
        state = tunnel_state(MACH, REYNOLDS_MAC, REFERENCE_MAC_M)
        for xc in TRIP_POSITIONS:
            for chord in (ROOT_CHORD_M, TIP_CHORD_M):
                values = []
                scale = math.sqrt(state["nu_inf_m2_s"] * xc * chord / state["U_inf_m_s"])
                for solution in (self.solution, dense):
                    profile = profile_metrics(solution["eta"] * scale,
                                              solution["u_over_U"] * state["U_inf_m_s"],
                                              state["nu_inf_m2_s"], state["U_inf_m_s"])
                    values.append(required_height(profile, 600., profile["y_m"][-1]))
                with self.subTest(xc=xc, chord=chord):
                    self.assertLess(abs(values[0] / values[1] - 1), 1e-5)

    def test_provided_reference_values_are_not_replaced_by_trapezoid_values(self):
        metadata = geometry_metadata()
        self.assertEqual(metadata["provided_reference_mac_m"], 4.629)
        self.assertEqual(metadata["provided_reference_area_m2"], 49.)
        self.assertAlmostEqual(metadata["linear_panel_area_m2"], 16.12776)
        self.assertAlmostEqual(metadata["linear_panel_trapezoid_mac_m"], 3.8082926829268287)
        for cfg, _, _ in self.results.values():
            self.assertEqual(cfg["wing"]["reference_length_m"], 4.629)
            self.assertEqual(cfg["wing"]["reference_area_m2"], 49.)
            self.assertEqual(cfg["conditions"][0]["reynolds"], 8e6)
            self.assertIn("physical AoA unspecified", cfg["conditions"][0]["alpha_label"])

    def test_all_positions_feasible_with_expected_uniform_heights(self):
        expected = {.07: .783788882, .10: .854548, .12: .893387}
        for xc, (_, summary, rows) in self.results.items():
            with self.subTest(xc=xc):
                self.assertEqual(summary["status"], "feasible")
                self.assertEqual(len(summary["design"]), 42)
                self.assertEqual(summary["uniform_height"]["status"], "feasible")
                self.assertAlmostEqual(summary["uniform_height"]["height_mm"], expected[xc], delta=2e-6)
                self.assertTrue(all(r["selected_re_k"] >= 600 for r in rows))
                # Upper and lower surfaces coincide only under the documented
                # identical flat-plate assumption, not as a real wing prediction.
                for eta in {r["eta"] for r in rows}:
                    pair = [r["selected_height_mm"] for r in rows if r["eta"] == eta]
                    self.assertEqual(pair[0], pair[1])

    def test_absolute_coordinates_follow_user_left_panel(self):
        for xc, (_, summary, _) in self.results.items():
            root = next(r for r in summary["design"] if r["eta"] == 0 and r["surface"] == "upper")
            tip = next(r for r in summary["design"] if r["eta"] == 1 and r["surface"] == "upper")
            self.assertEqual(root["y_span_m"], -1.333)
            self.assertEqual(tip["y_span_m"], -6.250)
            self.assertAlmostEqual(root["x_global_m"], 9.280 + xc * 5.56)
            self.assertAlmostEqual(tip["x_global_m"], 12.557 + xc)

    def test_generated_rows_and_line_endings(self):
        self.assertEqual(self.generated["profile_rows"], 21 * 3 * 2 * 1001)
        for filename in (*self.generated["config_files"], "geometry_metadata.json", "profiles.csv"):
            data = (self.directory / filename).read_bytes()
            self.assertNotIn(b"\r\n", data)
        metadata = json.loads((self.directory / "geometry_metadata.json").read_text(encoding="utf-8"))
        self.assertIn("not experimental confirmation", metadata["profile_laminar_confirmation"])


if __name__ == "__main__":
    unittest.main()
