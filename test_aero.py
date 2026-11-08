"""Numerical reference checks; run with python -m unittest test_aero -v."""

import unittest

import numpy as np

from aero import (build_surrogate_profile, edge_state, profile_metrics,
                  required_height, roughness_re, thwaites, tunnel_state)


class AerodynamicsTests(unittest.TestCase):
    def test_tunnel_re_and_cp_zero_recover_freestream(self):
        fs = tunnel_state(0.2, 3e6, .64607)
        self.assertAlmostEqual(fs["U_inf_m_s"] * .64607 / fs["nu_inf_m2_s"], 3e6)
        edge = edge_state(np.array([0.]), fs, .2)
        self.assertAlmostEqual(float(edge["U_e_m_s"][0]), fs["U_inf_m_s"], places=10)
        self.assertAlmostEqual(float(edge["nu_e_m2_s"][0]), fs["nu_inf_m2_s"], places=12)

    def test_cp_invalid_pressure_is_not_clipped(self):
        fs = tunnel_state(.2, 3e6, .6)
        for cp in [-100., 2.]:
            with self.assertRaises(ValueError):
                edge_state([cp], fs, .2)

    def test_thwaites_flat_plate_exact_integral(self):
        s = np.linspace(0, .5, 101)
        result = thwaites(s, np.full(s.shape, 50.), 1.5e-5)
        expected = np.sqrt(.45 * 1.5e-5 * s / 50.)
        np.testing.assert_allclose(result["theta_m"], expected, rtol=1e-12)
        # The .45 integral fit predicts .67082 rather than exact Blasius .664.
        self.assertLess(abs(result["theta_m"][-1] / (.664 * np.sqrt(1.5e-5 * .5 / 50)) - 1), .011)

    def test_downstream_start_requires_explicit_thickness(self):
        with self.assertRaises(ValueError):
            thwaites([.01, .02, .03], [20, 20, 20], 1.5e-5)
        theta0 = np.sqrt(.45 * 1.5e-5 * .01 / 20.)
        result = thwaites([.01, .02, .03], [20, 20, 20], 1.5e-5, theta0)
        self.assertAlmostEqual(result["theta_m"][-1], np.sqrt(.45 * 1.5e-5 * .03 / 20))

    def test_separation_history_is_retained(self):
        s = np.linspace(0, 1, 401)
        ue = 30 - 20 * np.minimum(s / .7, 1)
        result = thwaites(s, ue, 1.5e-5)
        self.assertIsNotNone(result["first_separation_s_m"])
        self.assertTrue(result["separated"][-1])
        self.assertTrue(np.isnan(result["H"][-1]))

    def test_quartic_matches_integral_targets(self):
        for lam in [-.06, 0, .08]:
            profile = build_surrogate_profile(1e-4, lam, 50, 1.5e-5)
            self.assertAlmostEqual(profile["theta_m"] / 1e-4, 1., places=5)
            self.assertAlmostEqual(profile["H"], profile["H_target"], places=5)
            self.assertLess(profile["delta99_m"], profile["y_m"][-1])

    def test_local_roughness_velocity_and_height(self):
        profile = profile_metrics([0, .001, .002], [0, 20, 20], 1e-5, 20)
        height = required_height(profile, 400, .001)
        self.assertAlmostEqual(height, np.sqrt(400 * 1e-5 / 20000), places=12)
        self.assertAlmostEqual(roughness_re(profile, height), 400., places=10)
        self.assertIsNone(required_height(profile, 400, .0001))
        with self.assertRaises(ValueError):
            roughness_re(profile, .003)
        with self.assertRaises(ValueError):
            required_height(profile, 400, .003)

    def test_variable_viscosity_uses_value_at_roughness_top(self):
        profile = profile_metrics([0, .001, .002], [0, 20, 20], [1e-5, 2e-5, 2e-5], 20)
        height = required_height(profile, 400, .001)
        self.assertGreater(height, np.sqrt(400 * 1e-5 / 20000))
        self.assertAlmostEqual(roughness_re(profile, height), 400., places=10)

    def test_bad_profiles_rejected(self):
        for y, u in [([.0001, .001, .002], [0, 20, 20]),
                     ([0, .001, .002], [1, 20, 20]),
                     ([0, .001, .002], [0, 20, 19]),
                     ([0, .001, .002], [0, 10, 15])]:
            with self.assertRaises(ValueError):
                profile_metrics(y, u, 1e-5, 20)


if __name__ == "__main__":
    unittest.main()
