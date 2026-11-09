"""Independent source-table and artifact checks for the exact six-series airfoil."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from generate_naca64a204 import (
    DEFAULT_OUTPUT, SOURCE_DATA, coordinates_and_provenance, generate,
    load_source_data,
)


class Naca64A204Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workspace = Path(__file__).resolve().parent
        cls.temporary = tempfile.TemporaryDirectory(prefix="naca64a204-tests-", dir=cls.workspace)
        cls.directory = Path(cls.temporary.name)
        cls.points, cls.provenance = coordinates_and_provenance()

    @classmethod
    def tearDownClass(cls):
        if cls.directory.resolve().parent != cls.workspace:
            raise RuntimeError("Refusing test cleanup outside the workspace")
        cls.temporary.cleanup()

    def test_all_coordinates_match_independent_published_table(self):
        # The table is independent of this Python translation and retains
        # the exact values printed by PDAS, before percent-to-chord conversion.
        self.assertLess(self.provenance["published_table_comparison_max_abs_chord"], 5.5e-7)
        self.assertEqual(load_source_data()["published_direct_table"][7],
                         [9.9473, 1.8746, 10.0527, -.6825])

    def test_shape_is_four_percent_64a_family_and_preserves_original_chord(self):
        self.assertEqual(self.points.shape, (401, 2))
        np.testing.assert_array_equal(self.points[[0, 200, -1]], [[1, 0], [0, 0], [1, 0]])
        self.assertAlmostEqual(self.provenance["computed_max_native_normal_thickness"],
                               .04, delta=2e-5)
        # The actual nose can be slightly ahead of the nominal camber-line LE.
        # Do not silently translate or rotate the source chord to that point.
        self.assertAlmostEqual(self.points[:, 0].min(), -6.637177625593017e-6, delta=1e-12)
        self.assertGreater(self.points[100, 1], 0)

    def test_reproduces_existing_sweep_airfoil_byte_for_byte(self):
        destination = self.directory / "generated"
        metadata = generate(destination)
        actual = (destination / "naca64a204.dat").read_bytes()
        expected = (DEFAULT_OUTPUT / "naca64a204.dat").read_bytes()
        self.assertEqual(actual, expected)
        self.assertNotIn(b"\r", actual)
        self.assertEqual(hashlib.sha256(actual).hexdigest(),
                         "1618a7a74cab75532e4402146691f5aef38ee3db5de3895de506e57eafa1682e")
        self.assertEqual(metadata["output_sha256"], hashlib.sha256(actual).hexdigest())
        self.assertEqual((destination / "source_data.json").read_bytes(), SOURCE_DATA.read_bytes())

    def test_generation_is_limited_to_source_artifacts_and_is_idempotent(self):
        destination = self.directory / "untouched_workflow"
        destination.mkdir()
        sentinel = b'{"user_workflow": true}\n'
        (destination / "workflow.json").write_bytes(sentinel)
        generate(destination)
        dat = destination / "naca64a204.dat"
        timestamp = dat.stat().st_mtime_ns
        generate(destination)
        self.assertEqual(dat.stat().st_mtime_ns, timestamp)
        self.assertEqual((destination / "workflow.json").read_bytes(), sentinel)
        self.assertEqual({p.name for p in destination.iterdir()},
                         {"workflow.json", "naca64a204.dat", "source_data.json", "airfoil_source.json"})
        metadata = json.loads((destination / "airfoil_source.json").read_text())
        self.assertIn("public domain", metadata["attribution"].lower())
        self.assertEqual(metadata["source_data_file"], "source_data.json")

    def test_corrupted_mapping_cannot_silently_generate_a_different_airfoil(self):
        source = copy.deepcopy(load_source_data())
        source["psi7"] = [v * 1.05 for v in source["psi7"]]
        with self.assertRaisesRegex(ValueError, "published-coordinate check failed"):
            coordinates_and_provenance(source)
        source = copy.deepcopy(load_source_data())
        source["eps7"].pop()
        malformed = self.directory / "malformed.json"
        malformed.write_text(json.dumps(source), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "201 finite"):
            load_source_data(malformed)


if __name__ == "__main__":
    unittest.main()
