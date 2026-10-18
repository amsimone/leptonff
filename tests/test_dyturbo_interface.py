"""Steering/contract tests. These do not pretend to validate DYTurbo physics."""
import copy
import io
import os
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'adapters/dyturbo'))
from run import MODES, steering, validate_result
from build import extract
from run_decay import decay_steering, ratio


class InterfaceTests(unittest.TestCase):
    def test_w_differential_steering(self):
        card=decay_steering(process='wp',calls=1000)
        for line in ('finord = false','intDimRes = 4','intDimCT = 6','intDimVJ = 7','makecuts = false','qed = false'):
            self.assertIn(line+'\n',card)
        with self.assertRaises(ValueError):decay_steering(process='z')
        with self.assertRaises(ValueError):decay_steering(calls=1001)

    def test_ratio_keeps_on_off_covariance(self):
        # Perfectly correlated numerator = 2*denominator has zero ratio error.
        mean=[3,6,0,0,0];cov=[0.]*25
        cov[0]=4;cov[1]=cov[5]=8;cov[6]=16
        self.assertEqual(ratio(mean,cov,1,0),{'value':2,'error':0})
        self.assertIsNone(ratio([0]*5,cov,1,0)['value'])

    def test_dependency_timestamps_and_internal_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            archive = base / 'dependency.tar.gz'
            with tarfile.open(archive, 'w:gz') as t:
                member = tarfile.TarInfo('dep/generated');member.size = 4;member.mtime = 12345
                t.addfile(member, io.BytesIO(b'data'))
                link = tarfile.TarInfo('dep/link');link.type = tarfile.SYMTYPE;link.linkname = 'generated'
                t.addfile(link)
            extract(archive, base / 'safe', 'dep')
            self.assertEqual((base / 'safe/dep/link').read_bytes(), b'data')
            self.assertFalse((base / 'safe/dep/link').is_symlink())
            self.assertEqual(os.stat(base / 'safe/dep/generated').st_mtime, 12345)
            with tarfile.open(archive, 'w:gz') as t:
                link.linkname = '../../outside'
                t.addfile(link)
            with self.assertRaisesRegex(RuntimeError, 'Escaping dependency'):
                extract(archive, base / 'unsafe', 'dep')

    def test_orders_and_methods(self):
        self.assertEqual(len(MODES), 6)
        for mode in MODES.values():
            card = steering(mode)
            for key in ('order_sudak', 'order_hcoef', 'order_ct', 'order_vjet'):
                self.assertIn(f'{key} = {mode.order}\n', card)
            self.assertIn('qed = false\n', card)
            self.assertIn('orderlha = false\n', card)
            self.assertIn(f'BORNquad = {str(mode.integration == "quadrature").lower()}\n', card)
        self.assertIn('fixedorder_only = false\n', steering(MODES['nlo-nll']))

    def test_fiducial_nnlo_guard(self):
        for name in ('nnlo-quad', 'nnlo-nnll'):
            with self.assertRaisesRegex(ValueError, 'cannot integrate NNLO'):
                steering(MODES[name], fiducial=True)
            card = steering(MODES[name], fiducial=True, hybrid=True)
            self.assertIn('VJquad = false\n', card)
            self.assertIn('BORNquad = true\n', card)
        steering(MODES['nnlo-vegas'], fiducial=True)
        steering(MODES['nlo-quad'], fiducial=True)

    def test_invalid_inputs(self):
        for kw in ({'mass': [0., 1.]}, {'qt': [1., 0.]}, {'calls': 0},
                   {'accuracy': float('nan')}, {'pdf': 'foo\nqed=true'}, {'process': 'h'}):
            with self.assertRaises(ValueError):
                steering(MODES['nlo-quad'], **kw)

    def test_mass_partition(self):
        card = steering(MODES['nlo-nll'], mass=(66., 68.), mass_bin_width=.5)
        self.assertIn('m_bins = [ 66 66.5 67 67.5 68 ]', card)
        with self.assertRaises(ValueError):
            steering(MODES['nlo-nll'], mass_bin_width=.7)

    def test_reject_silent_method_fallback_and_lost_order(self):
        mode = MODES['nlo-vegas']
        result = dict(schema='leptonff-dyturbo-native-1', units='fb', qcd_order=1,
                      fixed_order=True, fiducial=False,
                      resolved_orders=dict(hard=1, counterterm=1, vjet=1, sudakov=1),
                      bins=[dict(value=4., error=.5, terms=[
                          dict(name='signed', method='vegas', value=-2., error=.3),
                          dict(name='positive', method='vegas', value=6., error=.4)])])
        validate_result(result, mode)
        for mutate in (lambda r: r['resolved_orders'].update(vjet=0),
                       lambda r: r['bins'][0]['terms'][0].update(method='quadrature'),
                       lambda r: r['bins'][0].update(value=5.),
                       lambda r: r['bins'][0].update(error=.7)):
            changed = copy.deepcopy(result)
            mutate(changed)
            with self.assertRaises(ValueError):
                validate_result(changed, mode)


if __name__ == '__main__':
    unittest.main()
