"""End-to-end guards for the external laminar-integral sizing route."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from trip_design import integral_profile, run_design


class IntegralTripTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="integral_test_", dir=Path(__file__).resolve().parent)
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.row = dict(condition_id='a0', eta='0', surface='upper', x_c='0.1',
                        status='valid', reason='', theta_m='0.0001224744871391589',
                        shape_factor='2.6', nu_m2_s='0.000022685947564054', ue_m_s='68.057842692162',
                        edge_mach='0.2', laminar_confirmed='true', transition_xc='0.5', cf='0.002')
        self.config = dict(schema_version=1, name='Integral test',
                           wing=dict(semispan_m=1., sweep_le_deg=15., reference_length_m=1.,
                                     stations=[dict(eta=0., chord_m=1.)]),
                           conditions=[dict(id='a0', mach=.2, reynolds=3e6, alpha_deg=0.)],
                           trip=dict(type='grit', position_xc={'upper': .1}, height_margin=1.1,
                                     max_k_delta99=.8, available_heights_mm=[]),
                           analysis=dict(source='integral', integral_file='integral.csv'))

    def calculate(self, rows=None):
        rows = rows or [self.row]
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (self.directory / 'integral.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        config_path = self.directory / 'config.json'
        config_path.write_text(json.dumps(self.config))
        return run_design(config_path)

    def test_valid_integrals_remain_explicitly_a_surrogate(self):
        p = integral_profile([self.row], .1)
        self.assertIn('surrogate', p['method'])
        self.assertIsNone(p['lambda_thwaites'])
        self.assertAlmostEqual(p['theta_m']/float(self.row['theta_m']), 1., places=5)
        self.assertAlmostEqual(p['H'], 2.6, places=5)
        _, summary, conditions = self.calculate()
        self.assertEqual(summary['status'], 'feasible')
        self.assertGreaterEqual(conditions[0]['selected_re_k'], 600.)
        self.assertLessEqual(conditions[0]['selected_k_delta99'], .8)

    def test_transition_attachment_and_laminar_status_guards(self):
        for change in ({'transition_xc': '.1'}, {'transition_xc': '.09'},
                       {'transition_xc': 'NaN'}, {'cf': '0'}, {'cf': '-.001'},
                       {'laminar_confirmed': 'false'}, {'status': 'invalid', 'reason': 'solver failed'},
                       {'edge_mach': '.31'}, {'shape_factor': '4'}, {'theta_m': '-.001'}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                integral_profile([dict(self.row, **change)], .1)

    def test_exact_position_and_single_record_required(self):
        for rows in ([dict(self.row, x_c='.12')], [self.row, self.row]):
            with self.assertRaises(ValueError):
                integral_profile(rows, .1)

    def test_failed_condition_is_retained_and_invalidates_envelope(self):
        self.config['conditions'].append(dict(id='a2', mach=.2, reynolds=3e6, alpha_deg=2.))
        bad = dict(self.row, condition_id='a2', status='invalid', reason='no converged polar point',
                   theta_m='', shape_factor='', cf='', laminar_confirmed='false')
        _, summary, conditions = self.calculate([self.row, bad])
        self.assertEqual(len(conditions), 2)
        self.assertEqual(summary['design'][0]['status'], 'invalid')
        self.assertIn('no converged polar point', summary['design'][0]['reason'])
        self.assertIsNone(summary['design'][0]['height_mm'])

    def test_global_laminar_limit_applies_to_integral_route(self):
        self.config['analysis']['laminar_until_xc'] = .09
        _, summary, conditions = self.calculate()
        self.assertEqual(conditions[0]['status'], 'invalid')
        self.assertEqual(summary['status'], 'not_feasible')

    def test_high_freestream_mach_cannot_pass_through_low_edge_mach(self):
        self.config['conditions'][0]['mach'] = .5
        _, summary, conditions = self.calculate()
        self.assertEqual(conditions[0]['status'], 'invalid')
        self.assertIn('freestream Mach', conditions[0]['reason'])
        self.assertEqual(summary['status'], 'not_feasible')


if __name__ == '__main__':
    unittest.main()
