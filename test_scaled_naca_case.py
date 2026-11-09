"""Model dimensions and operating-point invariants for the requested 1:8 wing."""
import json
from pathlib import Path
import tempfile
import unittest

from run_scaled_naca64a204 import prepare
from section_workflow import discover, prepare_geometry


class ScaledNacaCaseTests(unittest.TestCase):
    def test_model_dimensions_and_reynolds_numbers_are_scaled_once(self):
        workspace = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(prefix='scaled-naca-test-', dir=workspace) as name:
            folder = Path(name)
            config = prepare(folder)
            settings = json.loads(config.read_text())
            self.assertEqual(settings['reference_mac_m'], .578625)
            self.assertEqual(settings['reference_area_m2'], .765625)
            self.assertEqual(len(settings['states']), 9)
            self.assertEqual({s['reynolds_mac'] for s in settings['states']}, {3e6, 5e6, 8e6})
            self.assertEqual({s['alpha_deg'] for s in settings['states']}, {0., 1., 2.})
            output = folder / 'output'
            output.mkdir()
            _, wing, geometry = prepare_geometry(settings, discover(folder / 'sections'), output)
            self.assertEqual(len(wing['stations']), 21)
            root, tip = wing['stations'][0], wing['stations'][-1]
            self.assertAlmostEqual(root['chord_m'], .695, places=10)
            self.assertAlmostEqual(tip['chord_m'], .125, places=10)
            self.assertEqual((root['y_m'], tip['y_m']), (-.166625, -.78125))
            self.assertAlmostEqual(root['x_le_m'], 1.16)
            self.assertAlmostEqual(tip['x_le_m'], 1.569625)
            self.assertAlmostEqual(wing['semispan_m'], .614625)
            self.assertAlmostEqual(8e6 * root['chord_m'] / settings['reference_mac_m'], 9608986.82220782)
            self.assertEqual(settings['alpha_reference'], 'global_x')
            self.assertAlmostEqual(geometry[0]['incidence_deg'], -.02233915, places=6)
            self.assertEqual(settings['xfoil']['panels'], 280)


if __name__ == '__main__':
    unittest.main()
