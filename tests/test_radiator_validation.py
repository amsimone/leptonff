"""Guards for the independent radiator and external normalized-order analysis."""
import importlib.util
import copy
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'adapters/validation'))
from check_mass_radiator import continuum, cumulative_coefficients, prediction
spec=importlib.util.spec_from_file_location('kkmc_analysis',ROOT/'adapters/kkmcee/compare.py')
kkmc=importlib.util.module_from_spec(spec);spec.loader.exec_module(kkmc)
sys.path.insert(0,str(ROOT/'adapters/kkmcee'))
from diagnose import fit_powers, validate_controls


class RadiatorTests(unittest.TestCase):
    def test_coupling_fit_recovers_terms_and_absolute_mc_covariance(self):
        factors=np.array([.0625,.125,.5,1,2])
        design=np.column_stack((factors,factors**2,factors**3))
        coefficients=np.array([0.,-60e-6,-27e-6])
        errors=np.array([5e-6,7e-6,5e-6,7e-6,11e-6])
        result=fit_powers(factors,design@coefficients,errors,[1,2,3])
        np.testing.assert_allclose(result['coefficients'],coefficients,atol=1e-17)
        expected=np.linalg.inv(design.T@np.diag(1/errors**2)@design)
        np.testing.assert_allclose(result['covariance'],expected,rtol=1e-11)
        # A perfect polynomial must not make the measured MC errors disappear.
        self.assertGreater(result['errors'][0],0)
        doubled=fit_powers(factors,design@coefficients,2*errors,[1,2,3])
        np.testing.assert_allclose(doubled['covariance'],4*expected,rtol=1e-11)

    def test_coupling_fit_rejects_unresolved_or_invalid_data(self):
        for errors in ([0.,1.],[1.,float('nan')]):
            with self.assertRaisesRegex(ValueError,'Invalid independent'):
                fit_powers([.5,1.],[1.,2.],errors,[2])
        with self.assertRaisesRegex(ValueError,'Unresolved'):
            fit_powers([1.,1.],[1.,1.],[1.,1.],[2,3])

    def test_residual_controls_require_independent_streams_and_fixed_settings(self):
        x=np.random.default_rng(42).normal(size=(20,10))
        central=self.sample_record(x,1)
        coupling=self.sample_record(x,2);coupling['alpha']/=2
        mass=self.sample_record(x,3);mass['mass']/=100
        validate_controls(central,[coupling],mass)
        mass['seed']=2
        with self.assertRaisesRegex(ValueError,'overlap'):validate_controls(central,[coupling],mass)
        mass['seed']=3;mass['cutoff']/=10
        with self.assertRaisesRegex(ValueError,'settings'):validate_controls(central,[coupling],mass)

    @staticmethod
    def sample_record(x,seed):
        return dict(schema='leptonff-kkmcee-fixed-q-v1',Q=91.1876,mass=.1056583755,
            alpha=1/137.035999084,cutoff=1e-5,fsr=True,modes=list(kkmc.MODES),
            trials=len(x),seed=seed,provenance=dict(client_sha256='client',defaults_sha256='defaults'),
            observables=[dict(observable='moment',value=2,mean=x.mean(axis=0).tolist(),
                covariance=(np.cov(x,rowvar=False)/len(x)).reshape(-1).tolist())])

    def test_unequal_seed_pool_matches_concatenated_trials(self):
        rng=np.random.default_rng(813)
        # Correlated weights, unequal chunk sizes and deliberately different
        # sample means exercise both within- and between-seed scatter.
        x=rng.normal(size=(73,10))+.1
        y=rng.normal(size=(127,10))+.2
        z=rng.normal(size=(59,10))+.05
        samples=[self.sample_record(a,s) for a,s in zip((x,y,z),(1,2,3))]
        result=kkmc.pool(samples)
        direct=self.sample_record(np.concatenate((x,y,z)),4)
        np.testing.assert_allclose(result['observables'][0]['mean'],direct['observables'][0]['mean'],atol=1e-15)
        np.testing.assert_allclose(result['observables'][0]['covariance'],direct['observables'][0]['covariance'],atol=1e-17)
        nested=kkmc.pool([kkmc.pool(samples[:2]),samples[2]])
        np.testing.assert_allclose(nested['observables'][0]['covariance'],direct['observables'][0]['covariance'],atol=1e-17)
        self.assertEqual(result['trials'],259)
        self.assertEqual(result['seeds'],[1,2,3])

    def test_pool_rejects_duplicate_seeds_and_incompatible_campaigns(self):
        x=np.random.default_rng(42).normal(size=(20,10))
        a=self.sample_record(x,1);b=self.sample_record(x,2)
        with self.assertRaisesRegex(ValueError,'overlapping seeds'):kkmc.pool([a,a])
        with self.assertRaisesRegex(ValueError,'overlapping seeds'):kkmc.pool([kkmc.pool([a,b]),a])
        for key,value in [('cutoff',1e-6),('alpha',.001),('mass',.2),('fsr',False)]:
            changed=copy.deepcopy(b);changed[key]=value
            with self.assertRaisesRegex(ValueError,'different physics'):kkmc.pool([a,changed])
        for key in ('client_sha256','defaults_sha256'):
            changed=copy.deepcopy(b);changed['provenance'][key]='other'
            with self.assertRaisesRegex(ValueError,'provenance'):kkmc.pool([a,changed])
        changed=copy.deepcopy(b);changed['observables'][0]['value']=3
        with self.assertRaisesRegex(ValueError,'different observables'):kkmc.pool([a,changed])
        changed=copy.deepcopy(b);changed['ceex_min_mass']=0
        with self.assertRaisesRegex(ValueError,'support floors'):kkmc.pool([a,changed])

    def test_pool_preserves_paired_order_cancellation(self):
        def weights(w):
            x=np.zeros((len(w),10))
            for k in range(5):
                x[:,k]=w*(1+.01*k)
                x[:,k+5]=x[:,k]*(.9+.001*k)
            return x
        a=self.sample_record(weights(np.array([.8,.9,1.1,1.2])),1)
        b=self.sample_record(weights(np.array([.7,.8,1.2,1.3,1.4])),2)
        result=kkmc.difference(kkmc.pool([a,b])['observables'][0],4,3)
        self.assertAlmostEqual(result['value'],.001,places=14)
        self.assertLess(result['mc_error'],1e-8)

    def test_first_mass_moment_independent_rational_integrals(self):
        c=continuum(2)
        self.assertAlmostEqual(c['P'],-4/3,places=14)
        self.assertAlmostEqual(c['C'],22/9,places=14)
        self.assertAlmostEqual(continuum(1)['K'],0,places=14)
        # The published raw sector is NOT silently made number conserving.
        self.assertAlmostEqual(c['G1'],-.719311605774509,places=9)

    def test_quadratic_matching_changes_only_the_finite_second_coefficient(self):
        rows=[dict(Q=91.1876,mass=m,scale=91.1876,alpha=1/137.035999084,N=3)
              for m in [.05282918775,.1056583755,.211316751]]
        differences=[prediction(r)['quadratic2']-prediction(r)['nll2'] for r in rows]
        self.assertLess(max(differences)-min(differences),1e-11)
        self.assertGreater(differences[0],0)

    def test_cumulative_double_plus_has_endpoint_enhancement(self):
        ordinary=cumulative_coefficients(.9);endpoint=cumulative_coefficients(.999)
        self.assertGreater(endpoint['PP'],ordinary['PP'])
        self.assertTrue(all(np.isfinite(v) for v in endpoint.values()))

    def test_normalized_order_difference_retains_joint_covariance(self):
        # A common production-weight fluctuation must cancel from ratios of
        # inclusive and measured weights, including between perturbative orders.
        sample=np.array([.8,.9,1.1,1.2])
        x=np.zeros((4,10))
        for k in range(5):
            x[:,k]=sample*(1+.01*k)
            x[:,k+5]=x[:,k]*(.9+.001*k)
        b=dict(mean=x.mean(axis=0).tolist(),
               covariance=(np.cov(x,rowvar=False)/len(x)).reshape(-1).tolist())
        r=kkmc.difference(b,4,3)
        self.assertAlmostEqual(r['value'],.001,places=14)
        self.assertLess(r['mc_error'],1e-8)
        b['mean'][4]=0
        with self.assertRaisesRegex(ValueError,'inclusive normalization'):kkmc.normalized(b,4)


if __name__=='__main__':unittest.main()
