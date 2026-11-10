"""Model dimensions and operating-point invariants for the requested 1:8 wing."""
import json
import contextlib
import hashlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from run_scaled_naca64a204 import CASE, main, prepare, read_airfoil
from section_workflow import discover, prepare_geometry


class ScaledNacaCaseTests(unittest.TestCase):
    def test_custom_airfoil_preserves_source_and_uses_its_name(self):
        root = Path(__file__).resolve().parent
        source = root / 'examples/airfoils/naca2412.dat'
        source_bytes = source.read_bytes()
        with tempfile.TemporaryDirectory(prefix='custom-foil-test-', dir=root) as name:
            folder = Path(name)
            cfg = json.loads(prepare(folder, source).read_text())
            metadata = json.loads((folder / 'model_geometry.json').read_text())
            self.assertEqual(cfg['airfoil_name'], 'NACA 2412')
            self.assertNotIn('64A204', cfg['name'] + cfg['data_origin'])
            self.assertNotIn('PDAS', cfg['data_origin'])
            self.assertEqual((folder / 'source_airfoil.dat').read_bytes(), source_bytes)
            self.assertEqual(source.read_bytes(), source_bytes)
            self.assertEqual(metadata['airfoil_sha256'], hashlib.sha256(source_bytes).hexdigest())
            self.assertEqual(metadata['airfoil_file'], 'source_airfoil.dat')
            self.assertEqual(len(list((folder / 'sections').iterdir())), 21)

    def test_custom_cli_preparation_from_another_folder_with_spaces(self):
        root = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(prefix='custom-cli-test-', dir=root) as name:
            folder = Path(name)
            source = folder / 'my airfoil.dat'
            source.write_bytes((root / 'examples/airfoils/naca2412.dat').read_bytes())
            process = subprocess.run([sys.executable, str(root / 'run_scaled_naca64a204.py'),
                '--airfoil', 'my airfoil.dat', '--airfoil-name', 'My custom section',
                '--output', 'my results', '--prepare-only'], cwd=folder,
                capture_output=True, text=True, timeout=30)
            self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            settings = json.loads((folder / 'my results/inputs/workflow.json').read_text())
            workflow = json.loads((folder / 'my results/all_conditions/workflow_summary.json').read_text())
            self.assertEqual(settings['airfoil_name'], 'My custom section')
            self.assertEqual(len(workflow['solver_plan']), 189)
            self.assertTrue(all(p['mach'] == .2 for p in workflow['solver_plan']))
            self.assertFalse((folder / 'my results/all_conditions/boundary_layer.csv').exists())

    def test_custom_airfoil_cannot_replace_builtin_inputs_or_results(self):
        root = Path(__file__).resolve().parent
        source = root / 'examples/airfoils/naca2412.dat'
        with self.assertRaisesRegex(ValueError, 'separate case'):
            prepare(CASE, source)
        with patch('run_scaled_naca64a204.execute') as execute, contextlib.redirect_stderr(io.StringIO()):
            code = main(['--airfoil', str(source), '--output', str(root / 'results/left_wing_naca64a204_scale8')])
        self.assertEqual(code, 2)
        execute.assert_not_called()

    def test_summarize_and_airfoil_override_cannot_mix(self):
        with patch('run_scaled_naca64a204.execute') as execute, contextlib.redirect_stderr(io.StringIO()):
            code = main(['--airfoil', 'my.dat', '--summarize-only'])
        self.assertEqual(code, 2)
        execute.assert_not_called()

    def test_bad_custom_file_fails_before_generating_sections(self):
        root = Path(__file__).resolve().parent
        text = (root / 'examples/airfoils/naca2412.dat').read_text()
        variants = [text.split('\n', 1)[1], text.replace('1.0000', '5.0000'),
                    text.replace('0.0013', 'nan', 1), 'Name\n1 0\n0 0\n1 0\n']
        with tempfile.TemporaryDirectory(prefix='bad-foil-test-', dir=root) as name:
            folder = Path(name)
            for i, contents in enumerate(variants):
                source = folder / f'bad{i}.dat'
                source.write_text(contents)
                with self.subTest(i=i), self.assertRaises(ValueError):
                    prepare(folder / f'case{i}', source)
                self.assertFalse((folder / f'case{i}/sections').exists())

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
