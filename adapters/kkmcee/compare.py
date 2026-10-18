#!/usr/bin/env python3
"""Compare normalized KKMC fixed-Q observables with the audited leptonff radiator."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'validation'))
from check_mass_radiator import continuum, prediction

MODES=('eex1','eex2','ceex0','ceex1','ceex2')


def validate(d):
    if d.get('schema')!='leptonff-kkmcee-fixed-q-v1' or d['modes']!=list(MODES):
        raise ValueError('Not the supported fixed-Q KKMC reference')
    if any(d.get(k) for k in ('isr','ifi','weak','pairs','running_alpha')):
        raise ValueError('Not the fixed-alpha photonic FSR benchmark')
    if d['trials']<2 or d['Q']<=10*d['mass'] or d['alpha']<=0:
        raise ValueError('Invalid reference inputs')
    if not 0<=d.get('ceex_min_mass',1.)<d['Q']:
        raise ValueError('Invalid CEEX support floor')
    for b in d['observables']:
        x=np.asarray(b['mean']);c=np.asarray(b['covariance']).reshape(10,10)
        if len(x)!=10 or not np.all(np.isfinite(x)) or not np.all(np.isfinite(c)):
            raise ValueError('Invalid moments/covariance')
        if not np.allclose(c,c.T,rtol=1e-10,atol=1e-15) or min(np.diag(c)) < -1e-14:
            raise ValueError('Invalid covariance symmetry/diagonal')
    return d


def ratio(b,k):
    x=np.asarray(b['mean']);c=np.asarray(b['covariance']).reshape(10,10)
    if x[k]<=5*math.sqrt(max(0.,c[k,k])):
        raise ValueError('Unresolved inclusive normalization')
    g=np.zeros(10);g[k+5]=1/x[k];g[k]=-x[k+5]/x[k]**2
    return float(x[k+5]/x[k]),g,c


def normalized(b,k):
    value,g,c=ratio(b,k)
    return dict(value=value,mc_error=math.sqrt(max(0.,float(g@c@g))))


def difference(b,k,j):
    v,g,c=ratio(b,k);w,h,_=ratio(b,j)
    return dict(value=v-w,mc_error=math.sqrt(max(0.,float((g-h)@c@(g-h)))))


def read(path):
    return validate(json.loads(path.read_text()))


def seeds(d):
    return set(d['seeds'] if 'seeds' in d else [d['seed']])


def pool(samples):
    """Pool independent raw samples, including within- and between-seed scatter.

    The stored covariance is that of the mean, not of individual trials.
    Ratios and paired-order differences must be formed only after pooling.
    """
    if not samples:raise ValueError('No samples to pool')
    first=samples[0];seen=set()
    keys=('Q','mass','alpha','cutoff','fsr','modes','isr','ifi','weak','pairs',
          'running_alpha','remove_soft_photons','normalization')
    for d in samples:
        validate(d)
        if any(d.get(k)!=first.get(k) for k in keys):
            raise ValueError('Cannot pool different physics or cutoff settings')
        if d.get('ceex_active',d['fsr'])!=first.get('ceex_active',first['fsr']):
            raise ValueError('Cannot pool different CEEX settings')
        # The v1 legacy driver retained the pinned 1 GeV default without
        # exporting it. New records make that support restriction explicit.
        if d.get('ceex_min_mass',1.)!=first.get('ceex_min_mass',1.):
            raise ValueError('Cannot pool different CEEX support floors')
        for key in ('client_sha256','defaults_sha256'):
            if d.get('provenance',{}).get(key)!=first.get('provenance',{}).get(key):
                raise ValueError('Cannot pool different executable/defaults provenance')
        if seen & seeds(d):raise ValueError('Duplicate or overlapping seeds')
        seen.update(seeds(d))
        if [(b['observable'],b['value']) for b in d['observables']]!=[
                (b['observable'],b['value']) for b in first['observables']]:
            raise ValueError('Cannot pool different observables')
    if len(samples)==1:return first
    n=sum(d['trials'] for d in samples)
    result=dict(first,trials=n,seeds=sorted(seen),observables=[])
    result.pop('seed',None)
    result['seconds']=sum(d.get('seconds',0) for d in samples)
    result['max_resolved_record_missing_fraction']=max(
        d.get('max_resolved_record_missing_fraction',0) for d in samples)
    result['provenance']={k:first.get('provenance',{}).get(k)
                          for k in ('client_sha256','defaults_sha256')}
    result['provenance']['combination']='raw trial-weighted means and pooled sample scatter before normalization'
    result['samples']=[dict(trials=d['trials'],seeds=sorted(seeds(d))) for d in samples]
    for index,b in enumerate(first['observables']):
        mean=sum(d['trials']*np.asarray(d['observables'][index]['mean']) for d in samples)/n
        scatter=np.zeros((10,10))
        for d in samples:
            count=d['trials'];row=d['observables'][index]
            offset=np.asarray(row['mean'])-mean
            scatter+=count*(count-1)*np.asarray(row['covariance']).reshape(10,10)
            scatter+=count*np.outer(offset,offset)
        result['observables'].append(dict(observable=b['observable'],value=b['value'],
            mean=mean.tolist(),covariance=(scatter/(n*(n-1))).reshape(-1).tolist()))
    return validate(result)


def input_paths(folders):
    if isinstance(folders,Path):folders=[folders]
    paths=[]
    for folder in folders:
        found=sorted(folder.glob('*/result.json'))
        if not found:raise ValueError('No completed KKMC samples in '+str(folder))
        paths.extend(found)
    return paths


def records(folders):
    grouped={}
    for path in input_paths(folders):
        d=read(path);grouped.setdefault(d['Q'],[]).append(d)
    return [pool(grouped[q]) for q in sorted(grouped)]


def seed_consistency(folders):
    """Individual seed versus all other seeds, which are independent of it.

    These pulls overlap and observables are correlated: no global chi-square.
    """
    grouped={};rows=[]
    for path in input_paths(folders):
        d=read(path);grouped.setdefault(d['Q'],[]).append(d)
    for Q,samples in sorted(grouped.items()):
        for index,d in enumerate(samples):
            if len(samples)<2:continue
            rest=pool(samples[:index]+samples[index+1:])
            for b,c in zip(d['observables'],rest['observables']):
                x=normalized(b,4);y=normalized(c,4)
                error=math.hypot(x['mc_error'],y['mc_error'])
                rows.append(dict(Q=Q,seeds=sorted(seeds(d)),observable=b['observable'],
                    value=b['value'],difference=x['value']-y['value'],mc_error=error,
                    pull=(x['value']-y['value'])/error if error else 0))
    return rows


def same_inputs(a,b):
    return all(a[k]==b[k] for k in ('Q','mass','alpha','fsr')) and a['modes']==b['modes']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--central',type=Path,nargs='+',required=True)
    p.add_argument('--cutoff-control',type=Path,nargs='+',required=True)
    p.add_argument('--born-control',type=Path,nargs='+',required=True)
    p.add_argument('--zero-coupling-control',type=Path,nargs='+',required=True)
    p.add_argument('--small-alpha',type=Path,nargs='+',required=True)
    p.add_argument('--smaller-alpha',type=Path,nargs='+',required=True)
    p.add_argument('--profile',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    central=sorted(records(a.central),key=lambda d:d['Q'])
    if len({d['Q'] for d in central})!=len(central) or not central:
        raise ValueError('Need distinct parent masses')
    rows=[]
    for d in central:
        if not d['fsr']:raise ValueError('Central calculation has FSR off')
        for b in d['observables']:
            models={name:normalized(b,k) for k,name in enumerate(MODES)}
            row=dict(Q=d['Q'],observable=b['observable'],value=b['value'],kkmc=models,
                     ceex_order_change=difference(b,4,3),ceex_eex_difference=difference(b,4,1))
            if b['observable']=='moment':
                ref=prediction(dict(Q=d['Q'],mass=d['mass'],scale=d['Q'],alpha=d['alpha'],N=b['value']))
                row['leptonff']={k:ref[k] for k in ('order_alpha','ll','nll','quadratic')}
                row['scale_nll']=[prediction(dict(Q=d['Q'],mass=d['mass'],scale=d['Q']*s,
                                     alpha=d['alpha'],N=b['value']))['nll'] for s in (.5,1,2)]
                row['ceex2_minus_nll']=dict(value=models['ceex2']['value']-ref['nll'],
                                           mc_error=models['ceex2']['mc_error'])
            rows.append(row)
    controls={}
    cutoff=[]
    for d in records(a.cutoff_control):
        reference=next((c for c in central if same_inputs(c,d)),None)
        if reference is None or seeds(reference)&seeds(d) or reference['cutoff']==d['cutoff']:
            raise ValueError('Cutoff control must use the same physics and independent seeds')
        for b,c in zip(d['observables'],reference['observables']):
            if (b['observable'],b['value'])!=(c['observable'],c['value']):
                raise ValueError('Cutoff control observables differ')
            for k,name in enumerate(MODES):
                x=normalized(b,k);y=normalized(c,k)
                error=math.hypot(x['mc_error'],y['mc_error'])
                cutoff.append(dict(Q=d['Q'],observable=b['observable'],value=b['value'],mode=name,
                                   difference=x['value']-y['value'],mc_error=error,
                                   pull=(x['value']-y['value'])/error if error else 0))
    controls['cutoff']=cutoff
    controls['born']=[]
    for d in records(a.born_control):
        if d['fsr'] or d.get('ceex_active',False):
            raise ValueError('KKMC supports the no-radiation control in EEX, not CEEX')
        for b in d['observables']:
            for k in (0,1):
                x=normalized(b,k)
                if abs(x['value']-1)>1e-12:raise ValueError('Failed EEX FSR-off closure')
                controls['born'].append(dict(observable=b['observable'],value=b['value'],mode=MODES[k],result=x))
    controls['near_zero_coupling']=[]
    for d in records(a.zero_coupling_control):
        if d['alpha']>1e-7 or not d['fsr']:
            raise ValueError('Not the near-zero-coupling CEEX control')
        for b in d['observables']:
            x=normalized(b,4)
            if abs(x['value']-1)>max(1e-5,5*x['mc_error']):
                raise ValueError('Failed near-zero-coupling CEEX check')
            controls['near_zero_coupling'].append(dict(observable=b['observable'],value=b['value'],result=x))
    small=records(a.small_alpha);smaller=records(a.smaller_alpha)
    if len(small)!=1 or len(smaller)!=1 or seeds(small[0])&seeds(smaller[0]):
        raise ValueError('Need two independent coupling controls')
    lo,half=small[0],smaller[0]
    if (not math.isclose(lo['alpha'],2*half['alpha'],rel_tol=1e-13) or
            any(lo[k]!=half[k] for k in ('Q','mass','cutoff','fsr'))):
        raise ValueError('Coupling controls do not differ only in alpha')
    controls['first_order_limit']=[]
    for b,c in zip(lo['observables'],half['observables']):
        if b['observable']!='moment':continue
        x=normalized(b,4);y=normalized(c,4)
        aa=lo['alpha']/(2*math.pi);ah=half['alpha']/(2*math.pi)
        value=2*(y['value']-1)/ah-(x['value']-1)/aa
        error=math.hypot(2*y['mc_error']/ah,x['mc_error']/aa)
        moment=continuum(b['value']);L=2*math.log(lo['Q']/lo['mass'])
        expected=2*(L*moment['P']+moment['C'])
        controls['first_order_limit'].append(dict(N=b['value'],value=value,expected=expected,
            mc_error=error,pull=(value-expected)/error,
            caveat='linear alpha extrapolation; finite-alpha bias is not a statistical error'))
    profile=[json.loads(p.read_text()) for p in sorted(a.profile.glob('*.json'))]
    profile.sort(key=lambda d:d['cells'])
    if len(profile)!=2:raise ValueError('Need two fixed-Q profile resolutions')
    for b in profile[1]['rows']:
        previous=next(c for c in profile[0]['rows'] if c['N']==b['N'])
        row=next(r for r in rows if r['Q']==b['Q'] and r['observable']=='cumulative' and r['value']==b['cut'])
        row['leptonff']={name:b['bin_'+key] for name,key in
                        [('order_alpha','order_alpha'),('ll','ll'),('nll','nll'),('quadratic','quadratic')]}
        row['numerical_nll_shift']=abs(b['bin_nll']-previous['bin_nll'])
        row['ceex2_minus_nll']=dict(value=row['kkmc']['ceex2']['value']-b['bin_nll'],
                                   mc_error=row['kkmc']['ceex2']['mc_error'])
    a.output.mkdir(parents=True,exist_ok=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,len(central),figsize=(12,4),layout='constrained',sharey=True)
    for ax,d in zip(np.atleast_1d(axes),central):
        rr=[r for r in rows if r['Q']==d['Q'] and r['observable']=='moment']
        N=np.array([r['value'] for r in rr]);base=np.array([r['leptonff']['nll'] for r in rr])
        scale=np.array([r['scale_nll'] for r in rr])
        ax.fill_between(N,1e4*(scale.min(axis=1)/base-1),1e4*(scale.max(axis=1)/base-1),
                        color='#777777',alpha=.2,label='FF scale envelope')
        ax.plot(N,1e4*(np.array([r['leptonff']['ll'] for r in rr])/base-1),'o-',color='#2563eb',label='FF matched LL')
        ax.plot(N,1e4*(np.array([r['leptonff']['quadratic'] for r in rr])/base-1),'s--',color='#9333ea',label='FF quadratic finite')
        ax.errorbar(N,[1e4*r['ceex2_minus_nll']['value']/b for r,b in zip(rr,base)],
                    yerr=[1e4*r['ceex2_minus_nll']['mc_error']/b for r,b in zip(rr,base)],
                    fmt='o',color='black',label='KKMC CEEX2 (MC error)')
        ax.axhline(0,color='#dc2626',lw=1,label='FF matched NLL')
        ax.set_title(f'Q = {d["Q"]:g} GeV\n{d["trials"]/1e6:g}M weighted trials',fontsize=11)
        ax.set_xlabel(r'Moment index $N$');ax.grid(alpha=.2)
    np.atleast_1d(axes)[0].set_ylabel(r'$(R_N/R_N^{\mathrm{FF,NLL}}-1)\times10^4$')
    np.atleast_1d(axes)[0].legend(fontsize=7)
    fig.suptitle('Normalized bare-muon fixed-Q radiator; fixed α; ISR/IFI/weak/pairs off',fontsize=11)
    for suffix in ('png','pdf','svg'):fig.savefig(a.output/('fixed-q-moments.'+suffix),dpi=180)
    plt.close(fig)
    rr=[r for r in rows if r['observable']=='cumulative' and 'leptonff' in r]
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    xx=np.array([-math.log10(1-r['value']) for r in rr]);base=np.array([r['leptonff']['nll'] for r in rr])
    ax.errorbar(xx,[1e4*r['ceex2_minus_nll']['value']/b for r,b in zip(rr,base)],
                yerr=[1e4*r['ceex2_minus_nll']['mc_error']/b for r,b in zip(rr,base)],
                fmt='o',color='black',label='KKMC CEEX2 − FF NLL (MC error)')
    numeric=1e4*np.array([r['numerical_nll_shift'] for r in rr])/base
    ax.fill_between(xx,-numeric,numeric,color='#777777',alpha=.25,label='FF grid refinement shift')
    ax.plot(xx,1e4*(np.array([r['leptonff']['quadratic'] for r in rr])/base-1),'s--',
            color='#9333ea',label='FF quadratic finite − FF NLL')
    ax.axhline(0,color='#dc2626',lw=1);ax.grid(alpha=.2)
    ax.set_xlabel(r'$-\log_{10}(1-z_{\min})$, cumulative $z\geq z_{\min}$')
    ax.set_ylabel(r'Change relative to FF NLL [$10^{-4}$]')
    profile_Q=profile[1]['rows'][0]['Q']
    profile_trials=next(d['trials'] for d in central if d['Q']==profile_Q)
    ax.set_title(f'Q = {profile_Q:g} GeV; {profile_trials/1e6:g}M trials; photon-only radiator')
    ax.legend(fontsize=8)
    for suffix in ('png','pdf','svg'):fig.savefig(a.output/('fixed-q-cumulatives.'+suffix),dpi=180)
    plt.close(fig)
    moment_rows=[r for r in rows if r['observable']=='moment']
    cumulative_rows=[r for r in rows if r['observable']=='cumulative' and 'ceex2_minus_nll' in r]
    report=dict(schema='leptonff-kkmcee-comparison-v1',rows=rows,controls=controls,
                central_samples=[dict(Q=d['Q'],trials=d['trials'],seeds=sorted(seeds(d)),
                                      ceex_min_mass=d.get('ceex_min_mass',1.)) for d in central],
                cutoff_samples=[dict(Q=d['Q'],trials=d['trials'],seeds=sorted(seeds(d))) for d in records(a.cutoff_control)],
                seed_consistency=seed_consistency(a.central),
                seed_consistency_warning='CEEX2 seed-versus-rest pulls; overlapping tests and correlated observables, not a global chi-square.',
                interpretation='Resolved residuals are comparison diagnostics, not a certified missing-order uncertainty.',
                order_warning='CEEX1/CEEX2 and EEX1/EEX2 label YFS residual orders, not leptonff LL/NLL.',
                normalization_warning='Each prediction is divided by its own inclusive rate; CEEX excludes masses below its support floor (legacy default 1 GeV). No absolute partial-width comparison.',
                physical_validation_status='resolved residuals; their separation into omitted orders and sector conventions is not established',
                max_absolute_moment_difference=max(abs(r['ceex2_minus_nll']['value']) for r in moment_rows),
                max_absolute_moment_pull=max(abs(r['ceex2_minus_nll']['value']/r['ceex2_minus_nll']['mc_error']) for r in moment_rows),
                max_cutoff_pull=max(abs(r['pull']) for r in cutoff),
                max_cumulative_pull=max(abs(r['ceex2_minus_nll']['value']/r['ceex2_minus_nll']['mc_error']) for r in cumulative_rows),
                max_relative_cumulative_difference=max(abs(r['ceex2_minus_nll']['value']/r['leptonff']['nll']) for r in cumulative_rows),
                input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                    for folders in (a.central,a.cutoff_control,a.born_control,a.zero_coupling_control,
                                   a.small_alpha,a.smaller_alpha) for p in input_paths(folders)})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('max_absolute_moment_difference','max_absolute_moment_pull',
                     'max_cutoff_pull','max_cumulative_pull','max_relative_cumulative_difference')},indent=2))


if __name__=='__main__':main()
