"""Curve parsing, contour geometry, and unambiguous surface-mapping tests."""

from __future__ import annotations

import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

from curve_io import load_section, map_curve_to_section, parse_curve_name, read_curve


class CurveIOTests(unittest.TestCase):
    def setUp(self):
        self.workspace = Path(__file__).resolve().parent
        self.temporary = tempfile.TemporaryDirectory(prefix="curve-io-tests-", dir=self.workspace)
        self.directory = Path(self.temporary.name).resolve()
        self.assertEqual(self.directory.parent, self.workspace)

    def tearDown(self):
        if self.directory.resolve().parent != self.workspace:
            raise RuntimeError("Refusing cleanup outside test workspace")
        self.temporary.cleanup()

    def write_curve(self, points, name="section_state1_station-1.250", extra=""):
        path = self.directory / name
        lines = ["*KEYWORD", "*DEFINE_CURVE_TITLE", "test curve", "$ header", "1 0 1 1 0 0 0", "$ color", "$ X Y"]
        lines += [f" {x:.15E}  {y:.15E}" for x, y in points]
        path.write_text("\n".join(lines) + "\n*END\n" + extra, encoding="utf-8")
        return path

    @staticmethod
    def contour():
        return np.array([[1., 0.], [.75, .1], [.25, .12], [0., 0.],
                         [.25, -.1], [.75, -.06], [1., 0.]])

    @staticmethod
    def ellipse_contour():
        t = np.linspace(0, 2 * np.pi, 81)
        xy = np.column_stack((.5 + .5 * np.cos(t), .1 * np.sin(t)))
        xy[[0, 40, -1], 1] = 0.
        return xy

    def test_read_curve_accepts_fortran_notation_and_comments(self):
        path = self.write_curve(self.contour())
        lines = path.read_text().splitlines()
        lines[7] = " 1.0D+00 0.0d+00 $ inline comment"
        lines.insert(8, "$ whole-line comment")
        lines.insert(9, "")
        path.write_text("\n".join(lines) + "\n# permitted final comment\n")
        np.testing.assert_allclose(read_curve(path), self.contour(), atol=1e-15)

    def test_malformed_numeric_rows_and_multiple_blocks_rejected(self):
        for replacement in ("1 2 3", "NaN 2", "1 Inf", "not numeric", "*DEFINE_CURVE"):
            path = self.write_curve(self.contour())
            lines = path.read_text().splitlines()
            lines[8] = replacement
            path.write_text("\n".join(lines))
            with self.subTest(row=replacement), self.assertRaises(ValueError):
                read_curve(path)
        path = self.write_curve(self.contour(), extra="*KEYWORD\n*DEFINE_CURVE\n1 2\n*END\n")
        with self.assertRaisesRegex(ValueError, "after \\*END"):
            read_curve(path)
        path = self.write_curve(self.contour())
        path.write_text(path.read_text().replace("*END\n", ""))
        with self.assertRaisesRegex(ValueError, "missing terminal"):
            read_curve(path)

    def test_filename_metadata_and_decimal_suffix(self):
        self.assertEqual(parse_curve_name("cp_dist_state3_station-5.000"),
                         {"kind": "cp", "state": 3, "station_y_m": -5.})
        self.assertEqual(parse_curve_name("section_state12_station+1.25D+01.dat"),
                         {"kind": "section", "state": 12, "station_y_m": 12.5})
        for name in ("section_state0_station-1.0", "section_statex_station1", "cp_state1_station1", "section_state1_stationNaN"):
            with self.subTest(name=name), self.assertRaises(ValueError):
                parse_curve_name(name)

    def test_rotated_true_chord_normalization_and_surface_arc_lengths(self):
        angle = math.radians(17.)
        rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
        physical = 2. * self.contour() @ rotation.T + [5., -.3]
        section = load_section(self.write_curve(physical))
        self.assertAlmostEqual(section["chord_m"], 2., places=12)
        self.assertAlmostEqual(section["streamwise_chord_m"], 2 * math.cos(angle), places=12)
        self.assertAlmostEqual(section["incidence_deg"], 17., places=11)
        self.assertAlmostEqual(section["rotation_deg"], -17., places=11)
        np.testing.assert_allclose(section["normalized_x"], self.contour()[:, 0], atol=2e-14)
        np.testing.assert_allclose(section["normalized_z"], self.contour()[:, 1], atol=2e-14)
        upper = section["surfaces"]["upper"]
        self.assertEqual(upper["contour_indices"].tolist(), [3, 2, 1, 0])
        self.assertAlmostEqual(upper["s_m"][1], 2 * math.hypot(.25, .12), places=12)
        self.assertTrue(np.all(np.diff(upper["s_m"]) > 0))

    def test_geometry_labels_surfaces_when_first_branch_is_lower(self):
        points = self.contour()[::-1]
        section = load_section(self.write_curve(points))
        self.assertEqual(section["surfaces"]["upper"]["contour_branch"], 1)
        self.assertEqual(section["surfaces"]["lower"]["contour_branch"], 0)
        cp = np.column_stack((points[:, 0], np.arange(len(points), dtype=float)))
        mapped = map_curve_to_section(cp, section)
        np.testing.assert_array_equal(mapped["upper"]["cp"], [3, 4, 5, 6])
        np.testing.assert_array_equal(mapped["lower"]["cp"], [3, 2, 1, 0])

    def test_closed_te_and_repeated_le_preserve_raw_contour(self):
        points = np.insert(self.contour(), 4, [0., 0.], axis=0)
        section = load_section(self.write_curve(points))
        self.assertEqual(len(section["x_m"]), 8)
        self.assertEqual(section["leading_edge_indices"].tolist(), [3, 4])
        self.assertEqual(len(section["surfaces"]["upper"]["x_c"]), 4)
        self.assertEqual(len(section["surfaces"]["lower"]["x_c"]), 4)

    def test_nonmonotone_duplicate_and_crossing_geometry_rejected(self):
        bad_contours = []
        folded = self.contour()
        folded[1, 0] = .2
        bad_contours.append(folded)
        bad_contours.append(np.insert(self.contour(), 2, self.contour()[1], axis=0))
        crossing = self.contour()
        crossing[2, 1] = -.2
        bad_contours.append(crossing)
        touching = self.contour()
        touching[2, 1] = touching[4, 1]
        bad_contours.append(touching)
        for points in bad_contours:
            with self.subTest(points=points.tolist()), self.assertRaises(ValueError):
                load_section(self.write_curve(points))

    def test_exact_81_point_mapping_keeps_distinct_surface_cp(self):
        points = self.ellipse_contour()
        section = load_section(self.write_curve(points))
        values = np.column_stack((points[:, 0], np.arange(81, dtype=float)))
        result = map_curve_to_section(values, section)
        self.assertEqual(len(result["upper"]["cp"]), 41)
        self.assertEqual(len(result["lower"]["cp"]), 41)
        self.assertEqual(result["upper"]["cp"][10], 30)
        self.assertEqual(result["lower"]["cp"][10], 50)
        self.assertAlmostEqual(result["upper"]["x_c"][10], result["lower"]["x_c"][10])

    def test_80_point_curve_without_le_preserves_positive_start_arc_length(self):
        points = self.ellipse_contour()
        section = load_section(self.write_curve(points))
        values = np.column_stack((points[:, 0], np.arange(81, dtype=float)))
        values = np.delete(values, 40, axis=0)
        result = map_curve_to_section(values, section)
        self.assertEqual(len(result["upper"]["cp"]), 40)
        self.assertEqual(len(result["lower"]["cp"]), 40)
        self.assertEqual(result["upper"]["cp"][0], 39)
        self.assertEqual(result["lower"]["cp"][0], 41)
        self.assertGreater(result["upper"]["s_m"][0], 0.)
        self.assertGreater(result["lower"]["s_m"][0], 0.)

    def test_other_grids_interpolate_separate_surface_geometry(self):
        section = load_section(self.write_curve(self.contour()))
        values = np.array([[1., -1.], [.5, -2.], [0., 1.], [.5, .4], [1., 0.]])
        result = map_curve_to_section(values, section)
        np.testing.assert_allclose(result["upper"]["cp"], [1., -2., -1.])
        np.testing.assert_allclose(result["lower"]["cp"], [1., .4, 0.])
        self.assertAlmostEqual(result["upper"]["z_global_m"][1], .11)
        self.assertAlmostEqual(result["lower"]["z_global_m"][1], -.08)
        self.assertNotEqual(result["upper"]["s_m"][1], result["lower"]["s_m"][1])
        with self.assertRaisesRegex(ValueError, "extrapolation"):
            map_curve_to_section(values + [0.1, 0.], section)

    def test_pressure_conversion_requires_explicit_consistent_references(self):
        points = self.contour()
        section = load_section(self.write_curve(points))
        cp = np.arange(len(points), dtype=float) - 3.
        pressures = np.column_stack((points[:, 0], 101325. + 250. * cp))
        mapped = map_curve_to_section(pressures, section, "pressure", 250., 101325.)
        np.testing.assert_allclose(mapped["upper"]["cp"], cp[[3, 2, 1, 0]])
        for q, reference in ((None, 101325.), (250., None), (0., 101325.), (-1., 101325.), (250., float("nan"))):
            with self.subTest(q=q, reference=reference), self.assertRaises(ValueError):
                map_curve_to_section(pressures, section, "pressure", q, reference)
        pressures[1, 1] = -100.
        with self.assertRaisesRegex(ValueError, "Absolute pressure"):
            map_curve_to_section(pressures, section, "pressure", 250., 101325.)
        gauge = np.column_stack((points[:, 0], 250. * cp))
        mapped_gauge = map_curve_to_section(gauge, section, "gauge_pressure", 250., 0.)
        np.testing.assert_allclose(mapped_gauge["upper"]["cp"], cp[[3, 2, 1, 0]])


if __name__ == "__main__":
    unittest.main()
