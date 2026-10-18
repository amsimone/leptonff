"""Statistical and physics-content guards for independent Born-QCD FSR."""
import importlib.util
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
try:
    import numpy as np
except ImportError:
    np=None


def module(name,file):
    spec=importlib.util.spec_from_file_location(name,ROOT/file)
    value=importlib.util.module_from_spec(spec)
    sys.path.insert(0,str((ROOT/file).parent))
    try:spec.loader.exec_module(value)
    finally:sys.path.pop(0)
    return value


builder=module('fsr_build','adapters/renesance/build.py')
if np is not None:
    analysis=module('fsr_compare','adapters/renesance/compare.py')
    accuracy=module('fsr_accuracy','adapters/renesance/compare_leptonff.py')


class SourceGuardTests(unittest.TestCase):
    def test_schema_defaults_do_not_consume_nested_properties(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp);folder=source/'Generator/schema';folder.mkdir(parents=True)
            for name in ('renesance.schema','model.schema','foam.schema'):
                (folder/name).write_text('schema: {\nproperties: {\nfirst: {type: number, default: 1.25}\nsecond: {type: string, default: "CT14nlo"}\n}}\n')
            self.assertEqual(builder.defaults(source),{'first':'1.25','second':'CT14nlo'})

    def test_integrated_fsr_is_not_an_iqed_switch(self):
        # A pinned upstream integrated NC formula contains all three radiation
        # classes even when iqed=4. Its FSR terms depend on the current at s.
        terms=[]
        for name,count in [('v0s',4),('v0spr',2),('a2sspr',2),('vas',4),('vaspr',1),('a1sspr',1)]:
            terms.extend('     &   +'+name+'*alpha\n' for _ in range(count))
        text='      subroutine nc_ha_1313_1616_1spr (s,spr,hard)\n        hard=\n'+''.join(terms)+'\n        hard = hard*conhc*cfprime\n      end\n'
        selected=builder.neutral_fsr_integrated(text,'1313_1616')
        self.assertEqual(selected.count('+v0s*'),4);self.assertEqual(selected.count('+vas*'),4)
        self.assertNotIn('+v0spr',selected);self.assertNotIn('+a2sspr',selected)
        with self.assertRaises(ValueError):builder.neutral_fsr_integrated(text.replace('+vas*alpha\n','',1),'1313_1616')


@unittest.skipIf(np is None,'Optional NumPy not installed')
class StatisticsTests(unittest.TestCase):
    def record(self,seed=1):
        d=dict(schema='leptonff-born-fsr-v1',born_qcd=True,pid=202,ecm=13000,alpha=.0073,
               muon_mass=.10566,parent_min=50,parent_max=200,pdf='common',pdf_member=0,
               weak=False,isr=False,ifi=False,pairs=False,running_alpha=False,
               mode='fo',integration='monte_carlo',trials=10,seed=seed,
               normalization='unitary at each parent Q; independent SANC inclusive counterterm',
               normalization_quadrature_error=1e-8,settings={k:'0' for k in analysis.PHYSICS},
               channels=['born','sv_delta','real','photos','pythia','our_order_alpha','our_ll','our_nll','inclusive_counterterm'])
        d['bins']=[dict(observable='fiducial',lower=0,upper=0,mean=[100,-20,10,0,0,0,0,0,0],covariance=[0]*81)]
        return d

    def test_correlated_on_off_ratio_cancels_source_noise(self):
        d=self.record();b=d['bins'][0]
        event=np.array([[80,-16,8],[90,-18,9],[110,-22,11],[120,-24,12]],dtype=float)
        x=np.zeros((4,9));x[:,:3]=event
        b['mean']=x.mean(axis=0).tolist();b['covariance']=(np.cov(x,rowvar=False)/len(x)).reshape(-1).tolist()
        r=analysis.ratio(b,{0:1,1:1,2:1})
        self.assertAlmostEqual(r['correction'],-.1,places=14)
        self.assertLess(r['error'],1e-8)
        self.assertGreater(r['born_error'],0)

    def test_independent_seeds_and_quadrature_cannot_be_faked(self):
        a=self.record(1);b=self.record(2)
        combined=analysis.combine([a,b]);self.assertEqual(combined['trials'],20)
        with self.assertRaises(ValueError):analysis.combine([a,a])
        a['integration']=b['integration']='quadrature'
        with self.assertRaises(ValueError):analysis.combine([a,b])

    def test_additive_matching_subtracts_the_entire_first_order(self):
        reference=self.record()['bins'][0]
        ff=self.record()['bins'][0]
        ff['mean']=[100,0,0,0,0,90,91,91.25,0]
        result=analysis.match_to_reference(ff,reference,7)
        self.assertAlmostEqual(result['correction'],-.0875,places=14)
        ff['mean'][7]=ff['mean'][5]
        self.assertAlmostEqual(analysis.match_to_reference(ff,reference,7)['correction'],
                               analysis.ratio(reference,{0:1,1:1,2:1})['correction'],places=14)
        cov=np.zeros((9,9));cov[5,5]=cov[7,7]=cov[5,7]=cov[7,5]=4
        ff['covariance']=cov.reshape(-1).tolist()
        # An all-mode run has the joint covariance: its identical O(alpha)
        # contributions must cancel as random variables, not just as means.
        self.assertLess(analysis.match_to_reference(ff,ff,7)['error'],1e-12)

    def test_mismatched_physics_and_numerical_controls_rejected(self):
        a=self.record();b=self.record(2);b['pairs']=True
        with self.assertRaises(ValueError):analysis.combine([a,b])
        b=self.record(2);b['parent_max']=400
        with self.assertRaises(ValueError):analysis.combine([a,b])
        b=self.record(2);b['settings']['ome']='1e-6'
        with self.assertRaises(ValueError):analysis.combine([a,b])

    def test_accuracy_difference_retains_shared_covariance(self):
        b=self.record()['bins'][0]
        x=np.zeros((4,9))
        x[:,0]=[80,90,110,120]
        x[:,5]=.90*x[:,0];x[:,6]=.91*x[:,0];x[:,7]=.9125*x[:,0]
        b['mean']=x.mean(axis=0).tolist()
        b['covariance']=(np.cov(x,rowvar=False)/len(x)).reshape(-1).tolist()
        r=accuracy.paired_ratio(b,7,6)
        self.assertAlmostEqual(r['value'],.9125/.91-1,places=14)
        self.assertLess(r['mc_error'],1e-8)

    def scans(self):
        import copy
        scans=[]
        for cells,mass,y,shift in [(2048,32,64,0),(1024,32,64,.001),(2048,16,32,.002)]:
            scan={}
            for factor in (.5,1.,2.):
                d=self.record()
                d.update(mode='analytic',integration='quadrature',ff_bins=cells,
                         frag_scale_factor=factor,frag_scale=91.1876*factor)
                d['settings'].update(ff_bins=str(cells),mass_nodes=str(mass),
                                     rapidity_nodes=str(y),frag_scale_factor=str(factor))
                d['bins'][0]['mean']=[100,0,0,0,0,90,91,91.25+shift,0]
                scan[factor]=copy.deepcopy(d)
            scans.append(scan)
        return scans

    def test_accuracy_controls_and_first_order_scale_cancellation(self):
        scans=self.scans()
        rows,checks=accuracy.analyse(scans)
        self.assertEqual(checks['max_order_alpha_scale_spread'],0)
        self.assertAlmostEqual(rows[0]['numerical']['nll_ll']['estimate'],.003/91,places=14)
        scans[0][.5]['bins'][0]['mean'][5]+=1e-5
        with self.assertRaisesRegex(ValueError,'first-order scale cancellation'):
            accuracy.analyse(scans)

    def test_accuracy_envelope_cannot_mix_resolutions_or_physics(self):
        scans=self.scans();scans[1][1.]['settings']['mass_nodes']='16'
        with self.assertRaisesRegex(ValueError,'separate coarser-grid'):
            accuracy.analyse(scans)
        scans=self.scans();scans[1][1.]['alpha']*=2
        with self.assertRaisesRegex(ValueError,'physical prediction'):
            accuracy.analyse(scans)


if __name__=='__main__':unittest.main()
