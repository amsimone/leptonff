#!/usr/bin/env python3
"""Leptonff-only Z mass comparison, separating QED scale and numerical effects."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from compare import COMMON, PHYSICS, validate

CHANNELS = (5, 6, 7)
NAMES = ('order_alpha', 'matched_ll', 'matched_nll')
METRICS = ('fo_born', 'll_born', 'nll_born', 'll_fo', 'nll_fo', 'nll_ll')
PAIRS = ((5, 0), (6, 0), (7, 0), (6, 5), (7, 5), (7, 6))
SCALES = (.5, 1., 2.)


def paired_ratio(b, numerator, denominator):
    """Fractional difference, with covariance of the two channels."""
    m = np.asarray(b['mean'], dtype=float)
    c = np.asarray(b['covariance'], dtype=float).reshape(len(m), len(m))
    den = m[denominator]
    if den <= 0 or den <= 5*math.sqrt(max(0., c[denominator, denominator])):
        raise ValueError('Unresolved/nonpositive ratio denominator')
    g = np.zeros(len(m))
    g[numerator] += 1/den
    g[denominator] -= m[numerator]/den**2
    return dict(value=float(m[numerator]/den-1),
                mc_error=math.sqrt(max(0., float(g@c@g))))


def signature(d):
    return ([d[k] for k in COMMON], {k:d['settings'][k] for k in PHYSICS},
            [(b['observable'], b['lower'], b['upper']) for b in d['bins']])


def load_scan(folder):
    """One deterministic calculation per scale; never average scale variations."""
    scan = {}
    for path in sorted(folder.glob('*.json')):
        d = json.loads(path.read_text())
        if d.get('schema') != 'leptonff-born-fsr-v1':
            continue
        validate(d)
        if d['pid'] != 202 or d['mode'] != 'analytic' or d['integration'] != 'quadrature':
            raise ValueError('This comparison requires deterministic, analytic Z mass records')
        factor = d.get('frag_scale_factor', 1.)
        if factor in scan:
            raise ValueError('Duplicate scale: quadrature repetitions are not statistics')
        if not math.isfinite(factor) or factor <= 0:
            raise ValueError('Invalid fragmentation scale')
        # Only explicitly scanned records may supply noncentral scales.
        if factor!=1. and ('frag_scale' not in d or 'frag_scale_factor' not in d['settings']):
            raise ValueError('Noncentral scale needs explicit run metadata')
        if 'frag_scale' in d and not math.isclose(d['frag_scale'], factor*float(d['settings']['mz']),
                                                rel_tol=1e-13):
            raise ValueError('Scale metadata does not match the declared pole-mass convention')
        if any(np.any(np.asarray(b['covariance'])) for b in d['bins']):
            raise ValueError('Deterministic records must not advertise MC covariance')
        d['_path'] = str(path)
        d['_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        scan[factor] = d
    if not scan:
        raise ValueError('No analytic Z records')
    central = scan.get(1.)
    if central is None or any(signature(d) != signature(central) for d in scan.values()):
        raise ValueError('Scale scan changes physics inputs or histogram definitions')
    numerical = ('ff_bins', 'mass_nodes', 'rapidity_nodes')
    if any(any(d['settings'].get(k) != central['settings'].get(k) for k in numerical)
           for d in scan.values()):
        raise ValueError('Scale envelope mixes numerical resolutions')
    fixed = {k:v for k,v in central['settings'].items() if k!='frag_scale_factor'}
    if any({k:v for k,v in d['settings'].items() if k!='frag_scale_factor'} != fixed
           for d in scan.values()):
        raise ValueError('Scale envelope changes other settings')
    return scan


def metrics(b):
    return {k:paired_ratio(b, *p) for k,p in zip(METRICS, PAIRS)}


def analyse(scans):
    """Require independent grid and quadrature refinements at every QED scale."""
    fine, grid, quad = scans
    if any(set(s) != set(SCALES) for s in scans):
        raise ValueError('Each scan must contain exactly scales 0.5, 1, 2')
    for factor in SCALES:
        f, g, q = (s[factor] for s in scans)
        if signature(f) != signature(g) or signature(f) != signature(q):
            raise ValueError('Numerical controls change the physical prediction')
        sf, sg, sq = (d['settings'] for d in (f,g,q))
        if not (int(sg['ff_bins']) < int(sf['ff_bins']) == int(sq['ff_bins']) and
                all(sf[k] == sg[k] for k in ('mass_nodes','rapidity_nodes')) and
                all(int(sq[k]) < int(sf[k]) for k in ('mass_nodes','rapidity_nodes'))):
            raise ValueError('Need separate coarser-grid and coarser-quadrature controls')
        # All remaining settings must agree (notably PDFs, widths and parent domain).
        allowed = {'ff_bins','mass_nodes','rapidity_nodes'}
        for control in (g,q):
            if {k:v for k,v in sf.items() if k not in allowed} != {
                    k:v for k,v in control['settings'].items() if k not in allowed}:
                raise ValueError('Numerical control changes other settings')
    nominal = fine[1.]
    rows = []
    oa_spread = born_spread = 0.
    for i,b in enumerate(nominal['bins']):
        if b['observable'] not in ('mll','fiducial'):
            continue
        fs = {s:metrics(fine[s]['bins'][i]) for s in SCALES}
        controls = {}
        for key in METRICS:
            grid_shift = max(abs(fs[s][key]['value']-metrics(grid[s]['bins'][i])[key]['value'])
                             for s in SCALES)
            quad_shift = max(abs(fs[s][key]['value']-metrics(quad[s]['bins'][i])[key]['value'])
                             for s in SCALES)
            controls[key] = dict(grid_shift=grid_shift, quadrature_shift=quad_shift,
                                 estimate=grid_shift+quad_shift)
        born = b['mean'][0]
        for s in SCALES:
            born_spread = max(born_spread, abs(fine[s]['bins'][i]['mean'][0]/born-1))
            oa_spread = max(oa_spread, abs(fs[s]['fo_born']['value']-fs[1.]['fo_born']['value']))
        rows.append(dict(observable=b['observable'], lower=b['lower'], upper=b['upper'],
                         rates_pb={name:b['mean'][ch] for name,ch in zip(('born',*NAMES),(0,*CHANNELS))},
                         central=fs[1.],
                         scale={key:dict(lower=min(fs[s][key]['value'] for s in SCALES),
                                         upper=max(fs[s][key]['value'] for s in SCALES),
                                         values={str(s):fs[s][key]['value'] for s in SCALES})
                                for key in METRICS}, numerical=controls))
    if born_spread > 1e-13 or oa_spread > 2e-12:
        raise ValueError(f'Failed first-order scale cancellation: {oa_spread:g}, Born {born_spread:g}')
    return rows, dict(max_order_alpha_scale_spread=oa_spread, max_born_scale_spread=born_spread)


def plot(rows, output, window=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bins = [r for r in rows if r['observable']=='mll' and
            (window is None or r['lower']>=window[0] and r['upper']<=window[1])]
    edges = np.array([r['lower'] for r in bins]+[bins[-1]['upper']])
    widths = np.diff(edges)
    fig, axes = plt.subplots(4,1,figsize=(9,10),sharex=True,layout='constrained',
                             gridspec_kw={'height_ratios':[2.1,1.2,1.3,1.3]})
    colors = ('#303030','#2563eb','#dc2626')
    labels = (r'$\mathcal{O}(\alpha)$',r'$\mathcal{O}(\alpha)+\mathrm{LL}$',
              r'$\mathcal{O}(\alpha)+\mathrm{NLL}_{\rm phot.}$')
    def band(ax, key, color, multiplier=100):
        y = np.array([r['central'][key]['value'] for r in bins])*multiplier
        lower = np.array([r['scale'][key]['lower'] for r in bins])*multiplier
        upper = np.array([r['scale'][key]['upper'] for r in bins])*multiplier
        ax.fill_between(edges,np.r_[lower,lower[-1]],np.r_[upper,upper[-1]],
                        step='post',color=color,alpha=.18,lw=0)
        ax.stairs(y,edges,color=color,lw=1.3,baseline=None)
        return y
    for color,label,name,key in zip(colors,labels,NAMES,('fo_born','ll_born','nll_born')):
        rate = np.array([r['rates_pb'][name] for r in bins])/widths
        axes[0].stairs(rate,edges,color=color,label=label,lw=1.4,baseline=None)
        lower = np.array([r['rates_pb']['born']*(1+r['scale'][key]['lower']) for r in bins])/widths
        upper = np.array([r['rates_pb']['born']*(1+r['scale'][key]['upper']) for r in bins])/widths
        axes[0].fill_between(edges,np.r_[lower,lower[-1]],np.r_[upper,upper[-1]],
                             step='post',color=color,alpha=.15,lw=0)
        band(axes[1],key,color)
    axes[0].stairs(np.array([r['rates_pb']['born'] for r in bins])/widths,
                   edges,color='#888888',ls='--',label='Born (FSR off)',lw=1,baseline=None)
    for key,color in (('ll_fo',colors[1]),('nll_fo',colors[2])):
        band(axes[2],key,color)
    y=band(axes[3],'nll_ll',colors[2],1000)
    num=np.array([r['numerical']['nll_ll']['estimate'] for r in bins])*1000
    axes[3].fill_between(edges,np.r_[y-num,(y-num)[-1]],np.r_[y+num,(y+num)[-1]],
                         step='post',facecolor='none',edgecolor='#555555',hatch='////',lw=.4,
                         label='Numerical refinement estimate')
    axes[0].set_yscale('log')
    axes[0].set_ylabel(r'$d\sigma/dm_{\mu\mu}$ [pb/GeV]')
    axes[1].set_ylabel(r'$(\sigma/\sigma_{\rm Born}-1)$ [%]')
    axes[2].set_ylabel(r'$(\sigma/\sigma_{\mathcal{O}(\alpha)}-1)$ [%]')
    axes[3].set_ylabel(r'$(\sigma_{\rm NLL}/\sigma_{\rm LL}-1)$ [‰]')
    axes[3].set_xlabel(r'Bare $m_{\mu\mu}$ [GeV]; bin width 0.25 GeV')
    axes[0].legend(fontsize=9,ncol=2)
    axes[3].legend(fontsize=8,loc='best')
    for ax in axes:
        ax.grid(alpha=.18)
    for ax in axes[1:]:
        ax.axhline(0,color='#777777',lw=.6)
    axes[0].set_title(r'$pp\to Z/\gamma^*\to\mu^+\mu^-$, 13 TeV; Born QCD; fixed $\alpha$, photons only'
                      '\nNo lepton cuts; parent Q = 50–200 GeV; unitary decay normalization',fontsize=10)
    fig.supxlabel(r'Shading: QED fragmentation-scale envelope, $\mu_{\rm frag}=(1/2,1,2)M_Z$.'
                  '\nDeterministic quadrature: no MC bands. Hatching: grid + quadrature shifts; no PDF/QCD scale variations.',
                  fontsize=9)
    stem = 'z-mll-leptonff'+('-peak' if window else '')
    for suffix in ('png','pdf','svg'):
        fig.savefig(output/(stem+'.'+suffix),dpi=180)
    plt.close(fig)
    return stem+'.png'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,required=True,help='Finest scale scan')
    p.add_argument('--grid-control',type=Path,required=True)
    p.add_argument('--quadrature-control',type=Path,required=True)
    p.add_argument('--previous-central',type=Path,
                   help='Optional existing central result for an explicit refinement comparison')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    scans=[load_scan(folder) for folder in (a.inputs,a.grid_control,a.quadrature_control)]
    rows,checks=analyse(scans)
    previous_info=None
    if a.previous_central:
        old=json.loads(a.previous_central.read_text());validate(old)
        if (old['integration']!='quadrature' or old['mode']!='analytic' or
                old.get('frag_scale_factor',1.)!=1. or signature(old)!=signature(scans[0][1.])):
            raise ValueError('Previous central record has incompatible physics')
        changes={key:0. for key in METRICS}
        for r,b in zip(rows,[b for b in old['bins'] if b['observable'] in ('fiducial','mll')]):
            previous=metrics(b)
            for key in METRICS:
                changes[key]=max(changes[key],abs(r['central'][key]['value']-previous[key]['value']))
        previous_info=dict(path=str(a.previous_central),
            sha256=hashlib.sha256(a.previous_central.read_bytes()).hexdigest(),max_ratio_shifts=changes)
    a.output.mkdir(parents=True,exist_ok=False)
    # Integrated peak ratios are evaluated from integrals, never averaged bin ratios.
    peak={}
    for s,d in scans[0].items():
        selected=[b for b in d['bins'] if b['observable']=='mll' and b['lower']>=90 and b['upper']<=92]
        mean=np.sum([b['mean'] for b in selected],axis=0)
        b=dict(mean=mean,covariance=[0]*81)
        key=str(float(s))
        peak[key]={name:float(mean[ch]) for name,ch in zip(('born',*NAMES),(0,*CHANNELS))}
        peak[key]['ratios']=metrics(b)
    report=dict(schema='leptonff-z-accuracy-comparison-v1',qcd='Born',scales=list(SCALES),
                convention='unitary, fixed-alpha photonic; pole-mass fragmentation scale',
                scale_uncertainty='three-point envelope, not a statistical confidence interval',
                numerical_uncertainty='sum of absolute grid and quadrature refinement shifts; not a rigorous bound',
                inputs=[{str(float(s)):dict(path=d['_path'],sha256=d['_sha256']) for s,d in scan.items()} for scan in scans],
                checks=checks,peak_90_92=peak,rows=rows)
    report['plots']=[plot(rows,a.output),plot(rows,a.output,(84,96))]
    report['checks']['max_nll_ll_numerical_estimate']=max(
        r['numerical']['nll_ll']['estimate'] for r in rows)
    report['checks']['resolved_nll_ll_bins']=sum(
        abs(r['central']['nll_ll']['value'])>10*r['numerical']['nll_ll']['estimate']
        for r in rows if r['observable']=='mll')
    report['checks']['mass_bins']=sum(r['observable']=='mll' for r in rows)
    if previous_info:
        report['previous_central']=previous_info
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    with (a.output/'z-mll-leptonff.tsv').open('w') as stream:
        writer=csv.writer(stream,delimiter='\t')
        writer.writerow(['observable','lower','upper','metric','central','MC_error',
                         'scale_min','scale_max','grid_shift','quadrature_shift'])
        for r in rows:
            for key in METRICS:
                writer.writerow([r['observable'],r['lower'],r['upper'],key,
                                 r['central'][key]['value'],r['central'][key]['mc_error'],
                                 r['scale'][key]['lower'],r['scale'][key]['upper'],
                                 r['numerical'][key]['grid_shift'],r['numerical'][key]['quadrature_shift']])
    print(json.dumps(dict(checks=report['checks'],peak_90_92=peak['1.0'],plots=report['plots']),indent=2))


if __name__=='__main__':
    main()
