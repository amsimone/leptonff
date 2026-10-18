"""Independent checks of conservative source reconstruction and region covariance."""
from pathlib import Path
import sys
import unittest
try:
    import numpy as np
except ImportError:
    np=None
if np is not None:
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'adapters/dyturbo'))
    from compare_precision import density_map,region_covariance


@unittest.skipIf(np is None,'Optional NumPy not installed')
class PrecisionTests(unittest.TestCase):
    def test_conservative_stencil_and_input_derivative(self):
        weights=np.array([1.,2.,4.,3.,5.,7.,11.]);h=.25
        mapping=density_map(weights,h)
        average=np.array([1,0,1/12,0,1/80])
        # Every reconstructed cell integral is exactly its own native bin.
        np.testing.assert_allclose(h*np.einsum('k,jki->ji',average,mapping),np.eye(len(weights)),atol=1e-13)
        c=mapping@weights
        shifted=weights.copy();shifted[2]+=.01
        np.testing.assert_allclose(mapping@shifted-c,.01*mapping[:,:,2],atol=1e-12)

    def test_region_recovers_negative_cross_bin_products(self):
        # A parent stratum starts in bin 1 and may migrate to bins 0/1/2.
        after=np.array([0,1,2,0,2,1]);n=len(after)
        basis=np.array([[t**k for k in range(5)] for t in np.linspace(-.4,.4,n)])
        events=np.zeros((n,3,5))
        for j,b in enumerate(after):events[j,b]+=basis[j];events[j,1]-=basis[j]
        mean=events.mean(axis=0)
        cov=np.array([np.cov(events[:,b,:],rowvar=False,ddof=1)/n for b in range(3)])
        for mask in (np.array([True,True,False]),np.array([True,False,True]),np.ones(3,dtype=bool)):
            recovered=region_covariance(mean,cov,n,mask,bool(mask[1]))
            direct=np.cov(events[:,mask,:].sum(axis=1),rowvar=False,ddof=1)/n
            np.testing.assert_allclose(recovered,direct,atol=1e-15)


if __name__=='__main__':unittest.main()
