#!/usr/bin/env python3
"""Keep numerical, scale, matching and external-comparison diagnostics separate."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'renesance'))
from compare_leptonff import load_scan, signature


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-scales',type=Path,required=True)
    p.add_argument('--variant-scales',type=Path,required=True)
    p.add_argument('--numerics',type=Path,required=True)
    p.add_argument('--external-report',type=Path,required=True)
    p.add_argument('--audit-report',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    baseline=load_scan(a.baseline_scales);variant=load_scan(a.variant_scales)
    if set(baseline)!={.5,1,2} or set(variant)!=set(baseline):
        raise ValueError('Need the same three scales for both prescriptions')
    for s in baseline:
        b,v=baseline[s],variant[s]
        if b.get('mass_matching','single_insertion')!='single_insertion' or v.get('mass_matching')!='quadratic_finite':
            raise ValueError('Unexpected matching prescriptions')
        if signature(b)!=signature(v):
            raise ValueError('Matching variation changes physics or bins')
        if {k:x for k,x in b['settings'].items() if k!='mass_matching'} != {
                k:x for k,x in v['settings'].items() if k!='mass_matching'}:
            raise ValueError('Matching variation changes other settings')
    numerics=json.loads(a.numerics.read_text());external=json.loads(a.external_report.read_text())
    audit=json.loads(a.audit_report.read_text())
    if not audit['passed']:raise ValueError('The independent coefficient audit failed')
    numerical={(r['observable'],r['lower'],r['upper']):r for r in numerics['rows']}
    rows=[];first_order=0.
    for i,b in enumerate(baseline[1]['bins']):
        if b['observable'] not in ('mll','fiducial'):continue
        nominal=b['mean'][7]
        values=[baseline[s]['bins'][i]['mean'][7]/nominal-1 for s in (.5,1,2)]
        completion=[variant[s]['bins'][i]['mean'][7]/nominal-1 for s in (.5,1,2)]
        for s in baseline:
            bb,vv=baseline[s]['bins'][i],variant[s]['bins'][i]
            first_order=max(first_order,abs((vv['mean'][5]-bb['mean'][5])/bb['mean'][0]))
        r=numerical[(b['observable'],b['lower'],b['upper'])]
        rows.append(dict(observable=b['observable'],lower=b['lower'],upper=b['upper'],
            nominal_nll_pb=nominal,scale_min=min(values),scale_max=max(values),
            quadratic_shift=completion[1],quadratic_scale_shifts=completion,
            diagnostic_envelope=[min(values+completion),max(values+completion)],
            numerical_estimate=r['numerical']['nll_born']['estimate']/(nominal/b['mean'][0])))
    if first_order>2e-12:raise ValueError('Matching variation changes the first-order result')
    selected=[i for i,b in enumerate(baseline[1]['bins']) if b['observable']=='mll' and
              b['lower']>=90 and b['upper']<=92]
    total=sum(baseline[1]['bins'][i]['mean'][7] for i in selected)
    peak=dict(nll_pb=total,
        scale_shifts={str(s):sum(baseline[s]['bins'][i]['mean'][7] for i in selected)/total-1 for s in baseline},
        quadratic_shift=sum(variant[1]['bins'][i]['mean'][7] for i in selected)/total-1)
    a.output.mkdir(parents=True,exist_ok=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for window in (None,(84,96)):
        rr=[r for r in rows if r['observable']=='mll' and
            (window is None or r['lower']>=window[0] and r['upper']<=window[1])]
        x=np.array([(r['lower']+r['upper'])/2 for r in rr])
        fig,ax=plt.subplots(figsize=(9,4.8),layout='constrained')
        ax.fill_between(x,[100*r['scale_min'] for r in rr],[100*r['scale_max'] for r in rr],
                        color='#777777',alpha=.25,label=r'QED fragmentation scale: $M_Z/2,M_Z,2M_Z$')
        ax.plot(x,[100*r['quadratic_shift'] for r in rr],color='#9333ea',
                label='Quadratic finite insertion / default NLL − 1')
        n=np.array([100*r['numerical_estimate'] for r in rr])
        ax.fill_between(x,-n,n,facecolor='none',edgecolor='#555555',hatch='////',lw=.4,
                        label='Default NLL numerical refinement estimate')
        ax.axhline(0,color='#dc2626',lw=1,label='Default matched NLL')
        ax.set_xlabel(r'Bare $m_{\mu\mu}$ [GeV]; 0.25 GeV bins')
        ax.set_ylabel('Change relative to default NLL [%]')
        ax.grid(alpha=.2);ax.legend(fontsize=8,loc='best')
        ax.set_title('Z mass: separate scale and matching diagnostics; Born QCD, photons only\n'
                     'Matching variation preserves NLL; no combined theory-error claim',fontsize=11)
        stem='z-mll-uncertainty-diagnostics'+('-peak' if window else '')
        for suffix in ('png','pdf','svg'):fig.savefig(a.output/(stem+'.'+suffix),dpi=180)
        plt.close(fig)
    report=dict(schema='leptonff-z-uncertainty-diagnostics-v1',rows=rows,peak_90_92=peak,
        max_first_order_matching_change=first_order,
        retained_coefficient_audit=dict(checks=len(audit['checks']),passed=audit['passed']),
        component_interpretation={
            'numerical':'empirical grid+quadrature refinement estimate, not a rigorous bound',
            'scale':'three-point QED fragmentation-scale envelope, not a confidence interval',
            'matching':'one explicitly NNLL finite-insertion completion; does not span all missing terms',
            'external':'KKMC fixed-Q residuals are observable-dependent comparisons, not mll error bands',
            'normalization':'conditional photon-only shape; an absolute partial-width correction is not supplied',
            'excluded_physics':'pairs, fermion-loop running, ISR, IFI, weak corrections and improved QCD production',
            'total':'no certified total theoretical uncertainty is assigned'},
        external_summary={k:external[k] for k in
            ('max_absolute_moment_difference','max_absolute_moment_pull','max_cutoff_pull',
             'max_cumulative_pull','max_relative_cumulative_difference','physical_validation_status')},
        input_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                      for path in (a.numerics,a.external_report,a.audit_report)})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(peak_90_92=peak,coefficient_audit=report['retained_coefficient_audit'],
                         external=report['external_summary']),indent=2))


if __name__=='__main__':main()
