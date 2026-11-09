"""XFOIL protocol/parser safeguards; tests do not require an installed solver."""
from __future__ import annotations

from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from xfoil_backend import (parse_cp, parse_dump, parse_polar, run_xfoil,
                           sample_laminar_station, xfoil_commands)


POLAR = """XFOIL Version 6.99
 Mach = 0.200 Re = 1.000 e 6 Ncrit = 9.000
 alpha CL CD CDp CM Top_Xtr Bot_Xtr
 ------ -------- --------- --------- -------- -------- --------
 2.000 .2000 .00600 .00100 -.0100 .6000 .7000
"""
X = np.array([1., .6, .2, .1, .01, .01, .1, .2, .6, 1.])
Y = np.array([0., .04, .06, .04, .001, -.001, -.04, -.06, -.04, 0.])
UE = np.array([1., 1.1, 1.2, 1.2, .1, -.1, -1.1, -1.1, -1., -.9])


def dump_text(kinematic=False, wake=True):
    rows = ["# s x y Ue/Vinf Dstar Theta Cf H"]
    for i, (x, y, ue) in enumerate(zip(X, Y, UE)):
        msq = .2**2 * ue**2 / (1 + .2 * .2**2 * (1 - ue**2))
        shape = (2.6 - .29 * msq) / (1 + .113 * msq) if kinematic else 2.6
        rows.append(f"{i*.2:.5f} {x:.5f} {y:.5f} {ue:.5f} .002600 .001000 .003000 {shape:.4f} 1.2 3 4 5")
    if wake:
        rows.append("1.80000 1.00010 .00000 1.00000 .002600 .001000 .000000 2.6000")
    return "\n".join(rows) + "\n"


class XfoilBackendTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="xfoil_test_", dir=Path(__file__).resolve().parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.airfoil = self.root / "input with spaces.dat"
        # Independent normalized closed contour; no numerical solver run.
        self.airfoil.write_text("_test\n" + "".join(f"{x} {y}\n" for x, y in zip(
            [1,.8,.6,.4,.2,0,.2,.4,.6,.8,1], [0,.03,.05,.06,.04,0,-.04,-.06,-.05,-.03,0])))
        self.executable = self.root / "xfoil fake.exe"
        self.executable.write_bytes(b"test fixture; never executed")

    def file(self, name, text):
        path = self.root / name
        path.write_text(text)
        return path

    def write_outputs(self, directory, *, polar=POLAR):
        (directory / "polar.txt").write_text(polar)
        (directory / "cp.txt").write_text("# x Cp\n" + "".join(f"{x} -.2\n" for x in X))
        (directory / "bl.txt").write_text(dump_text())

    def successful_process(self, args, **kwargs):
        self.write_outputs(Path(kwargs["cwd"]))
        self.assertFalse(kwargs["shell"])
        self.assertIn("XTR 1 1", kwargs["input"])
        return subprocess.CompletedProcess(args, 0, "viscous calculation complete")

    def run_mock(self, callback=None, **kwargs):
        with patch("xfoil_backend.subprocess.run", side_effect=callback or self.successful_process):
            return run_xfoil(self.executable, self.airfoil, self.root / "runs", 1e6, .2, 2., **kwargs)

    def test_cp_two_and_three_column_versions(self):
        for width in (2, 3):
            path = self.file(f"cp{width}.txt", "# x y Cp\n" + "".join(
                f"{x} {y} -.2\n" if width == 3 else f"{x} -.2\n" for x,y in zip(X,Y)))
            result = parse_cp(path)
            np.testing.assert_allclose(result["cp"], -.2)
            self.assertEqual(result["y_c"] is None, width == 2)

    def test_dump_mass_and_kinematic_h_formats(self):
        for kin in (False, True):
            # Wake row repeats rounded TE s, and has fewer columns than surface.
            result = parse_dump(self.file(f"dump{kin}.txt", dump_text(kin)), .2)
            self.assertEqual(len(result["x_c"]), 11)
            msq = .2**2 * UE**2 / (1 + .2 * .2**2 * (1 - UE**2))
            np.testing.assert_allclose(result["h_kinematic"][:10],
                                       (2.6 - .29*msq)/(1+.113*msq), atol=5.1e-5)
            self.assertIn("kinematic" if kin else "mass", result["h_column_interpretation"])

    def test_malformed_nonfinite_and_missing_rows_rejected(self):
        for content in ("# x Cp\n.1 nan\n", "# x Cp\n.1 -.1\n.2 *****\n", "# x Cp\n"):
            with self.assertRaises(ValueError):
                parse_cp(self.file("bad_cp.txt", content))
        with self.assertRaises(ValueError):
            parse_dump(self.file("bad_bl.txt", dump_text().replace(".001000", "nan", 1)), .2)

    def test_polar_requires_one_matching_converged_point_and_conditions(self):
        result = parse_polar(self.file("polar.txt", POLAR), 2)
        self.assertEqual(result["top_xtr"], .6)
        self.assertEqual(result["reynolds"], 1e6)
        for content in (POLAR.rsplit("\n", 2)[0]+"\n", POLAR + "2 .2 .006 .001 -.01 .6 .7\n",
                        POLAR.replace("2.000", "3.000"), POLAR.replace(".6000", "NaN")):
            with self.assertRaises(ValueError):
                parse_polar(self.file("bad_polar.txt", content), 2)

    def test_success_logs_and_fresh_directories(self):
        first, second = self.run_mock(), self.run_mock()
        self.assertTrue(first["converged"], first["reason"])
        self.assertNotEqual(first["workdir"], second["workdir"])
        self.assertTrue(Path(first["stdout_path"]).is_file())
        self.assertTrue(Path(first["commands_path"]).is_file())
        self.assertIn("XTR 1 1", Path(first["commands_path"]).read_text())

    def test_timeout_nonzero_and_unconverged_are_failures(self):
        def timeout(*args, **kwargs):
            raise subprocess.TimeoutExpired(args[0], 1, output=b"partial log")
        result = self.run_mock(timeout)
        self.assertFalse(result["converged"])
        self.assertIn("timed out", result["reason"])
        self.assertEqual(Path(result["stdout_path"]).read_text(), "partial log")
        for code, output in [(1,"crashed"),(0,"VISCAL: Convergence failed")]:
            result = self.run_mock(lambda args, **kw: subprocess.CompletedProcess(args, code, output))
            self.assertFalse(result["converged"])

    def test_cp_dump_alone_do_not_establish_convergence(self):
        def no_polar(args, **kwargs):
            self.write_outputs(Path(kwargs["cwd"]), polar=POLAR.rsplit("\n", 2)[0]+"\n")
            return subprocess.CompletedProcess(args, 0, "")
        self.assertFalse(self.run_mock(no_polar)["converged"])

    def test_recovered_march_warning_requires_final_converged_polar(self):
        warning = 'MRCHDU: Convergence failed at 6 side 1\n'
        def recovered(args, **kwargs):
            self.write_outputs(Path(kwargs['cwd']))
            return subprocess.CompletedProcess(args, 0, warning + 'Point added to stored polar')
        result = self.run_mock(recovered)
        self.assertTrue(result['converged'], result['reason'])
        self.assertEqual(result['solver_warnings'], [warning.strip()])
        def no_final_point(args, **kwargs):
            self.write_outputs(Path(kwargs['cwd']), polar=POLAR.rsplit('\n', 2)[0]+'\n')
            return subprocess.CompletedProcess(args, 0, warning)
        self.assertFalse(self.run_mock(no_final_point)['converged'])

    def test_numerical_retries_stop_at_first_converged_solution(self):
        from unittest.mock import patch
        from xfoil_backend import run_xfoil_with_retries
        failures = [dict(converged=False, reason='not converged', panels=280, iterations=200),
                    dict(converged=True, reason='', panels=280, iterations=600)]
        with patch('xfoil_backend.run_xfoil', side_effect=failures) as solve:
            result = run_xfoil_with_retries('exe', 'foil', 'output', 5e6, .2, 2.,
                                           panels=280, retry_panels=[280, 240, 200])
        self.assertTrue(result['converged'])
        self.assertEqual(len(result['attempts']), 2)
        self.assertEqual(solve.call_count, 2)
        self.assertEqual(solve.call_args_list[0].args, solve.call_args_list[1].args)
        self.assertEqual(solve.call_args_list[1].kwargs['iterations'], 600)
        with patch('xfoil_backend.run_xfoil', return_value=failures[1]) as solve:
            run_xfoil_with_retries('exe', 'foil', 'output', 5e6, .2, 2., retry_panels=[240])
        self.assertEqual(solve.call_count, 1)

    def test_wrong_polar_conditions_rejected(self):
        def wrong_re(args, **kwargs):
            self.write_outputs(Path(kwargs["cwd"]), polar=POLAR.replace("1.000 e 6", "2.000 e 6"))
            return subprocess.CompletedProcess(args, 0, "")
        self.assertFalse(self.run_mock(wrong_re)["converged"])

    def test_laminar_sample_and_duplicate_printed_le_x(self):
        result = self.run_mock()
        for side in ("upper", "lower"):
            point = sample_laminar_station(result, side, .15)
            self.assertAlmostEqual(point["theta_c"], .001)
            self.assertGreater(point["ue_over_uinf"], 0)
            self.assertTrue(point["laminar_confirmed_by_model"])

    def test_transition_bracket_separation_and_high_local_mach_rejected(self):
        result = self.run_mock()
        with self.assertRaisesRegex(ValueError, "transition"):
            sample_laminar_station(result, "upper", .6)
        result["polar"]["top_xtr"] = .18
        with self.assertRaisesRegex(ValueError, "bracket"):
            sample_laminar_station(result, "upper", .15)
        result["polar"]["top_xtr"] = .6
        result["boundary_layer"]["cf"][4] = -.001
        with self.assertRaisesRegex(ValueError, "separated"):
            sample_laminar_station(result, "upper", .15)
        result["boundary_layer"]["cf"][4] = .003
        result["boundary_layer"]["ue_over_uinf_signed"][:5] *= 2
        with self.assertRaisesRegex(ValueError, "Mach"):
            sample_laminar_station(result, "upper", .15)

    def test_invalid_inputs_do_not_start_solver(self):
        with patch("xfoil_backend.subprocess.run") as process:
            for values in ((-1,.2,0), (1e6,.8,0), (1e6,.2,float('nan'))):
                with self.assertRaises(ValueError):
                    run_xfoil(self.executable, self.airfoil, self.root, *values)
            process.assert_not_called()

    def test_protocol_leaves_natural_transition_active(self):
        commands = xfoil_commands(1e6,.2,2)
        self.assertIn("PLOP\nG\n\n", commands)
        self.assertIn("VPAR\nN 9\nXTR 1 1\n", commands)
        self.assertIn("PACC\npolar.txt\n\nALFA 2\n", commands)


if __name__ == "__main__":
    unittest.main()
