#!/usr/bin/env python3
"""Compare FSR corrections on Born QCD, retaining paired ratio covariance."""
import argparse
import json
import math
from pathlib import Path
import numpy as np

COMMON = ('pid','ecm','alpha','muon_mass','parent_min','parent_max','pdf','pdf_member')
PHYSICS = ('mw','ww','mz','wz','igmz','irun','iscale','gfscheme')
CURVES = {'fo': ('SANC fixed O(α)', {0:1,1:1,2:1}),
          'photos': ('PHOTOS', {3:1}), 'pythia': ('PYTHIA (decay MEC)', {4:1}),
          'our_order_alpha': ('leptonff O(α)', {5:1}),
          'our_ll': ('leptonff matched LL', {6:1}),
          'our_nll': ('leptonff matched NLL', {7:1})}


def validate(data):
    if data.get('schema')!='leptonff-born-fsr-v1' or not data.get('born_qcd'):
        raise ValueError('Expected a Born-QCD FSR reference')
    if any(data.get(x) for x in ('weak','isr','ifi','pairs','running_alpha')):
        raise ValueError('Physics content differs from pure photonic FSR')
    n=len(data['channels'])
    if n!=9 or not data.get('normalization','').startswith('unitary at each parent Q'):
        raise ValueError('Missing independent inclusive normalization counterterm')
    if data['normalization_quadrature_error']>2e-6:
        raise ValueError('Inclusive normalization quadrature is unresolved')
    for b in data['bins']:
        x=np.asarray(b['mean']);c=np.asarray(b['covariance']).reshape(n,n)
        if x.shape!=(n,) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(c)):
            raise ValueError('Invalid bin values/covariance')
        scale=max(1,float(np.max(np.abs(c))))
        if not np.allclose(c,c.T,rtol=1e-10,atol=1e-12*scale) or np.min(np.diag(c)) < -1e-12*scale:
            raise ValueError('Invalid covariance symmetry/diagonal')
    return data


def combine(records):
    if not records:
        raise ValueError('No independent reference records')
    base=records[0]
    for d in records:
        validate(d)
        if any(d[k]!=base[k] for k in COMMON) or any(d['settings'][k]!=base['settings'][k] for k in PHYSICS):
            raise ValueError('Born physics inputs differ')
        if d['mode']!=base['mode'] or d['integration']!=base['integration']:
            raise ValueError('Cannot combine different algorithms')
        if d['settings']!=base['settings']:
            raise ValueError('Numerical settings differ: keep refinement controls separate')
        if [(b['observable'],b['lower'],b['upper']) for b in d['bins']]!=[(b['observable'],b['lower'],b['upper']) for b in base['bins']]:
            raise ValueError('Histogram definitions differ')
    if len({d['seed'] for d in records})!=len(records):
        raise ValueError('Duplicate random seeds')
    if base['integration']=='quadrature' and len(records)!=1:
        raise ValueError('Repeated quadrature is not independent statistics')
    weights=np.array([d['trials'] for d in records],dtype=float);weights/=weights.sum()
    result=dict(base);result['seeds']=[d['seed'] for d in records]
    result['trials']=sum(d['trials'] for d in records);result['bins']=[]
    for i,b in enumerate(base['bins']):
        mean=sum(w*np.asarray(d['bins'][i]['mean']) for w,d in zip(weights,records))
        cov=sum(w*w*np.asarray(d['bins'][i]['covariance']) for w,d in zip(weights,records))
        result['bins'].append(dict(b,mean=mean.tolist(),covariance=cov.tolist()))
    return result


def ratio(bin, coefficients):
    mean=np.asarray(bin['mean']);cov=np.asarray(bin['covariance']).reshape(len(mean),len(mean))
    B=mean[0]
    if B<=0 or B<=5*math.sqrt(max(0,cov[0,0])):
        return None
    coeff=np.zeros(len(mean))
    for k,v in coefficients.items():coeff[k]=v
    numerator=coeff@mean
    gradient=coeff/B;gradient[0]-=numerator/B**2
    var=float(gradient@cov@gradient)
    return {'correction':float(numerator/B-1),'error':math.sqrt(max(0,var)),
            'born':float(B),'born_error':math.sqrt(max(0,cov[0,0]))}


def match_to_reference(ff_bin, reference_bin, channel):
    """Additive matching to independent massive NLO; no noisy C1 repeated.

    This is a separately labelled prediction, not a replacement for the
    independent check of the library's own leading-power O(alpha) match.
    """
    if ff_bin is reference_bin:
        return ratio(ff_bin,{0:1,1:1,2:1,channel:1,5:-1})
    ref=ratio(reference_bin,{0:1,1:1,2:1})
    delta=ratio(ff_bin,{channel:1,5:-1})
    if ref is None or delta is None:return None
    return dict(ref,correction=ref['correction']+delta['correction']+1,
                error=math.hypot(ref['error'],delta['error']))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs',type=Path,nargs='+')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--match-to-reference',action='store_true',
                   help='Separate additive matching: SANC NLO + (FF resummed - FF NLO)')
    a=p.parse_args()
    paths=[]
    for path in a.inputs:paths.extend(sorted(path.glob('*.json')) if path.is_dir() else [path])
    grouped={}
    for path in paths:
        d=json.loads(path.read_text())
        if d.get('schema')!='leptonff-born-fsr-v1':continue
        validate(d);grouped.setdefault((d['pid'],d['mode']),[]).append(d)
    if not grouped:raise ValueError('No reference results')
    combined={k:combine(v) for k,v in grouped.items()}
    a.output.mkdir(parents=True,exist_ok=False)
    report={'schema':'leptonff-born-fsr-comparison-v1','processes':{},
            'matching':'SANC NLO + (FF resummed - FF O(alpha))' if a.match_to_reference else 'original library matching'}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=['#222222','#d97706','#7c3aed','#0891b2','#2563eb','#dc2626']
    for pid in sorted({k[0] for k in combined}):
        modes={k[1]:v for k,v in combined.items() if k[0]==pid}
        if a.match_to_reference and 'fo' in modes and 'analytic' in modes:
            reference,ff=modes['fo'],modes['analytic']
            if reference['integration']=='monte_carlo' and ff['integration']=='monte_carlo' and set(reference['seeds'])&set(ff['seeds']):
                raise ValueError('Independent additive matching needs disjoint MC seeds; use all mode for paired covariance')
        base=next(iter(modes.values()))
        for d in modes.values():
            if any(d[k]!=base[k] for k in COMMON) or any(d['settings'][k]!=base['settings'][k] for k in PHYSICS):
                raise ValueError('Comparison inputs differ')
        curves={}
        for name,(label,coeff) in CURVES.items():
            if a.match_to_reference and name=='our_order_alpha':continue
            mode='fo' if name=='fo' else 'showers' if name in ('photos','pythia') else 'analytic'
            if mode not in modes:mode='all'
            if mode not in modes:continue
            if name=='pythia' and not modes[mode]['pythia_mec']:label='PYTHIA (MEC off)'
            if a.match_to_reference and name in ('our_ll','our_nll'):
                label='SANC NLO + FF '+('LL' if name=='our_ll' else 'NLL')
            curves[name]=(modes[mode],label,coeff)
        def curve_ratio(name,d,index,coeff):
            if a.match_to_reference and name in ('our_ll','our_nll'):
                ref=modes.get('fo',modes.get('all'))
                if ref is None:raise ValueError('Additive matching requires an independent fixed-order result')
                return match_to_reference(d['bins'][index],ref['bins'][index],6 if name=='our_ll' else 7)
            return ratio(d['bins'][index],coeff)
        process='z' if pid==202 else 'wplus' if pid==222 else 'wminus'
        pr={'inputs':{k:base[k] for k in COMMON},'modes':{k:{'trials':v['trials'],'seeds':v['seeds'],
             'integration':v['integration']} for k,v in modes.items()},'rates':{},'plots':[]}
        for name,(d,label,coeff) in curves.items():
            pr['rates'][name]={b['observable']:curve_ratio(name,d,i,coeff) for i,b in enumerate(d['bins']) if b['observable'] in ('inclusive','fiducial')}
        for observable in ('mll',) if pid==202 else ('ptl','absetal','mt'):
            indices=[i for i,b in enumerate(base['bins']) if b['observable']==observable]
            if not indices:continue
            fig,axes=plt.subplots(2,1,figsize=(8,7),sharex=True,gridspec_kw={'height_ratios':[2,1.3]},layout='constrained')
            for color,(name,(d,label,coeff)) in zip(colors,curves.items()):
                bins=[d['bins'][i] for i in indices];centers=np.array([(b['lower']+b['upper'])/2 for b in bins]);widths=np.array([b['upper']-b['lower'] for b in bins])
                c=np.zeros(len(d['channels']))
                for k,v in coeff.items():c[k]=v
                rr=[curve_ratio(name,d,i,coeff) for i in indices]
                rates=np.array([(1+r['correction'])*r['born'] if r else np.nan for r in rr])/widths
                axes[0].plot(centers,rates,color=color,label=label,lw=1.3)
                delta=np.array([100*r['correction'] if r else np.nan for r in rr]);error=np.array([100*r['error'] if r else np.nan for r in rr])
                axes[1].plot(centers,delta,color=color,lw=1.2)
                if np.any(error>0):axes[1].fill_between(centers,delta-error,delta+error,color=color,alpha=.12,lw=0)
            centers=[(base['bins'][i]['lower']+base['bins'][i]['upper'])/2 for i in indices]
            axes[0].plot(centers,[base['bins'][i]['mean'][0]/(base['bins'][i]['upper']-base['bins'][i]['lower']) for i in indices],color='#888888',ls='--',label='Born, FSR off',lw=1)
            axes[0].set_ylabel('dσ/dO [pb / unit O]');axes[1].set_ylabel('FSR correction [%]');axes[1].axhline(0,color='gray',lw=.6)
            axes[1].set_xlabel({'mll':'bare m(μμ) [GeV]','ptl':'bare pT(μ) [GeV]','absetal':'bare |η(μ)|','mt':'mT(μ, ν) [GeV]'}[observable])
            for ax in axes:ax.grid(alpha=.2)
            axes[0].legend(fontsize=8,ncol=2)
            axes[0].set_title(f'{process}: Born QCD, pp 13 TeV; fixed α, photons only\n'+('no lepton cuts; parent Q = 50–200 GeV' if pid==202 else 'pT(μ) > 25 GeV, |η(μ)| < 2.4; parent Q = 60–200 GeV'),fontsize=10)
            if a.match_to_reference:fig.suptitle('Additive matching to independent massive NLO FSR',fontsize=10)
            stem=f'{process}-{observable}'
            for suffix in ('png','pdf','svg'):fig.savefig(a.output/(stem+'.'+suffix),dpi=180)
            plt.close(fig);pr['plots'].append(stem+'.png')
        with (a.output/(process+'-spectrum.tsv')).open('w') as stream:
            stream.write('observable\tlower\tupper\tcurve\tfsr_correction\tMC_error\n')
            for name,(d,_,coeff) in curves.items():
                for i,b in enumerate(d['bins']):
                    r=curve_ratio(name,d,i,coeff)
                    if r is not None:stream.write(f"{b['observable']}\t{b['lower']}\t{b['upper']}\t{name}\t{r['correction']:.16g}\t{r['error']:.16g}\n")
        report['processes'][process]=pr
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v['rates'] for k,v in report['processes'].items()},indent=2))


if __name__=='__main__':main()
