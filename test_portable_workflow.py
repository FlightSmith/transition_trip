"""Verify user-input overrides and launching from a different directory."""
import contextlib
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from generate_demo import generate
from run_trip import menu
from trip_design import main, run_design


class PortableWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root = Path(__file__).resolve().parent
        cls.temp = tempfile.TemporaryDirectory(prefix='portable-tests-', dir=cls.root)
        cls.folder = Path(cls.temp.name).resolve()
        assert cls.folder.parent == cls.root
        cls.inputs = cls.folder / 'inputs with spaces'
        generate(cls.inputs)

    @classmethod
    def tearDownClass(cls):
        if cls.folder.resolve().parent != cls.root:
            raise RuntimeError('Test cleanup escaped the workspace')
        cls.temp.cleanup()

    def test_cp_override_switches_profile_mode_and_preserves_configuration(self):
        path = self.inputs / 'config_profiles.json'
        original = path.read_bytes()
        cfg, summary, rows = run_design(path, cp_file=self.inputs / 'cp.csv', position_xc=.12)
        self.assertEqual(path.read_bytes(), original)
        self.assertEqual(cfg['analysis']['source'], 'cp')
        self.assertEqual(summary['source'], 'cp')
        self.assertEqual(summary['overrides'], {'source': 'cp', 'input_filename': 'cp.csv', 'position_xc': .12})
        self.assertTrue(all(r['x_c'] == .12 for r in rows))
        self.assertEqual(summary['inputs'][1]['sha256'], hashlib.sha256((self.inputs / 'cp.csv').read_bytes()).hexdigest())
        self.assertNotIn('supplied-profile', cfg['name'])
        self.assertIn('12% chord', cfg['name'])
        self.assertIn('template provenance does not establish', cfg['data_origin'])

    def test_profiles_override_switches_cp_mode(self):
        _, summary, rows = run_design(self.inputs / 'config.json', profiles_file=self.inputs / 'profiles.csv')
        self.assertEqual(summary['source'], 'profiles')
        self.assertEqual(summary['status'], 'feasible')
        self.assertTrue(all(r['profile_method'] == 'supplied laminar profile' for r in rows))

    def test_position_label_and_placeholder_aoa_survive_overrides(self):
        config = json.loads((self.inputs / 'config.json').read_text(encoding='utf-8'))
        config['conditions'][0]['alpha_label'] = 'placeholder; physical AoA unspecified'
        path = self.inputs / 'placeholder_config.json'
        path.write_text(json.dumps(config), encoding='utf-8')
        cfg, summary, rows = run_design(path, cp_file=self.inputs / 'cp.csv', position_xc=.12)
        self.assertIn('12% chord', cfg['name'])
        self.assertNotIn('10% chord', cfg['name'])
        self.assertTrue(any('physical AoA unspecified' in warning for warning in summary['warnings']))
        self.assertTrue(all('physical AoA unspecified' in r['alpha_label'] for r in rows
                            if r['condition_id'] == config['conditions'][0]['id']))

    def test_overrides_reject_invalid_position_and_conflicting_sources(self):
        path = self.inputs / 'config.json'
        for value in (0, 1, -.1, float('nan'), float('inf')):
            with self.subTest(value=value), self.assertRaises(ValueError):
                run_design(path, position_xc=value)
        with self.assertRaisesRegex(ValueError, 'either cp_file or profiles_file'):
            run_design(path, cp_file='a.csv', profiles_file='b.csv')

    def test_cli_works_from_other_directory_with_spaces(self):
        output = self.folder / 'results with spaces'
        process = subprocess.run([
            sys.executable, str(self.root / 'run_trip.py'),
            '--config', str(self.inputs / 'config_profiles.json'),
            '--cp-file', str(self.inputs / 'cp.csv'), '--position-xc', '.1',
            '--output', str(output), '--no-plot'], cwd=self.folder,
            capture_output=True, text=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
        self.assertTrue((output / 'report.html').is_file())
        summary = json.loads((output / 'summary.json').read_text(encoding='utf-8'))
        self.assertEqual(summary['source'], 'cp')
        self.assertEqual(summary['uniform_height']['height_mm'], .25)

    def test_missing_override_reports_input_error(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = main(['--config', str(self.inputs / 'config.json'),
                         '--cp-file', str(self.inputs / 'missing.csv'), '--no-plot',
                         '--output', str(self.folder / 'missing-result')])
        self.assertEqual(code, 2)
        self.assertFalse((self.folder / 'missing-result' / 'report.html').exists())

    def test_menu_custom_cp_preserves_quoted_paths(self):
        with patch('builtins.input', side_effect=['4', f'"{self.inputs / "cp.csv"}"', '.12',
                                                 f'"{self.folder / "custom result"}"']), contextlib.redirect_stdout(io.StringIO()):
            jobs = menu()
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0][jobs[0].index('--cp-file') + 1], str(self.inputs / 'cp.csv'))
        self.assertEqual(jobs[0][jobs[0].index('--position-xc') + 1], '0.12')


if __name__ == '__main__':
    unittest.main()
