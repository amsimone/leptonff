#!/usr/bin/env python3
"""Order normalization and independent-piece guards need no external generator."""
import copy
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'adapters/mcmule'))
from analyse import combine, normalized_coefficients, validate
from prepare_source import NEEDLE, REPLACEMENT
from summarize import fit


def record(seed=1):
    d = dict(schema='leptonff-mcmule-fixed-order-v1', Q=91.1876, mass=.1056583755,
             electron_mass=.00051099895069, support_min_mass=0, physical_alpha=1/137.035999084,
             internal_alpha=1, born_alpha_one_MeV_minus2=5e-10, xi=.2, mu_factor=1,
             softcut=1e-10, collcut=1e-14, current='pure vector', fsr=True, isr=False,
             ifi=False, pairs=False, weak=False, coupling='fixed', angular_map='fsr-direct',
             seed=seed, bad_phase_space_points=0, build_manifest_sha256='same',
             piece='ee2mmRRMMMM', measurement='shape',
             observables=[dict(observable='moment', value=n, integral=0. if n == 1 else -2.,
                               mc_error=0. if n == 1 else .1) for n in range(1, 15)])
    return d


class AnalysisTests(unittest.TestCase):
    def test_expanded_normalized_shape(self):
        # sigma=1+alpha*3+alpha²*7; M=1+alpha*5+alpha²*11.
        c = normalized_coefficients(2, 0, 4, 0, 3, 0)
        self.assertEqual(c['c1'], 2)
        self.assertEqual(c['c2'], -2)
        for alpha in (1e-4, 1e-5):
            exact = (1+5*alpha+11*alpha**2)/(1+3*alpha+7*alpha**2)
            expansion = 1+c['c1']*alpha+c['c2']*alpha**2
            self.assertLess(abs(exact-expansion), 10*alpha**3)

    def test_independent_product_variance(self):
        c = normalized_coefficients(2, .04, 4, .09, 3, .01)
        self.assertAlmostEqual(c['variance2'], .09+9*.04+4*.01+.01*.04)

    def test_seed_average(self):
        a, b = record(), record(2)
        b['observables'][1].update(integral=-4, mc_error=.2)
        validate([a, b])
        c = combine([a, b], 'ee2mmRRMMMM', 'shape', 1)
        self.assertEqual(c['value'], -3)
        self.assertAlmostEqual(c['variance'], (.1**2+.2**2)/4)
        self.assertEqual(c['seed_scatter_variance'], 1)

    def test_repeat_stream_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Reused'):
            validate([record(), record()])

    def test_mismatched_settings(self):
        for key, value in [('xi', .05), ('mass', .01), ('support_min_mass', 1),
                           ('angular_map', 'native'), ('physical_alpha', .008),
                           ('build_manifest_sha256', 'other')]:
            with self.subTest(key=key):
                a, b = record(), record(2); b[key] = value
                with self.assertRaises(ValueError):
                    validate([a, b])

    def test_invalid_piece_result(self):
        for key, value in [('internal_alpha', .007), ('bad_phase_space_points', 1),
                           ('measurement', 'unknown'), ('pairs', True),
                           ('running_alpha', True), ('piece', 'ee2mmAFMM')]:
            a = record(); a[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate([a])

    def test_invalid_measurements(self):
        a = record(); a['observables'][0]['integral'] = .1
        with self.assertRaisesRegex(ValueError, 'vanish'):
            validate([a])
        a = record(); a['measurement'] = 'inclusive'
        with self.assertRaisesRegex(ValueError, 'Inclusive'):
            validate([a])
        a = record(); a['observables'][1]['mc_error'] = math.nan
        with self.assertRaises(ValueError):
            validate([a])

    def test_dispatch_workaround_preserves_native_formula(self):
        self.assertEqual(REPLACEMENT.count('call psx2('), 1)
        self.assertIn('associated(gen, psx2)', REPLACEMENT)
        self.assertIn(NEEDLE.splitlines()[1].strip(), REPLACEMENT)
        self.assertNotIn('weight =', REPLACEMENT)

    def test_logarithmic_fit_and_absolute_covariance(self):
        L = [9., 11., 12., 13., 16., 18.]
        center = 2*math.log(91.1876/.1056583755)
        values = [-28+.5*(x-center)+.2*(x-center)**2 for x in L]
        result = fit(L, values, [1.]*6, [0, 1, 2])
        for actual, expected in zip(result['coefficients'], [-28, .5, .2]):
            self.assertAlmostEqual(actual, expected, places=11)
        self.assertEqual(result['dof'], 3)
        self.assertGreater(result['errors'][0], 0)
        scaled = fit(L, values, [2.]*6, [0, 1, 2])
        self.assertAlmostEqual(scaled['covariance'][0][0], 4*result['covariance'][0][0])

    def test_near_to_soft_settings_cannot_be_pooled(self):
        a, b = record(), record(2)
        a['nts_switch'] = .001;b['nts_switch'] = .0001
        with self.assertRaisesRegex(ValueError, 'near-to-soft'):
            validate([a, b])


if __name__ == '__main__':
    unittest.main()
