#!/usr/bin/env python3
"""Recover common-weight rate estimates from the optional raw iteration audit.

MCFM's per-histogram inverse-variance result is reconstructed independently.
Equal weights are applied within each fixed-call QCD contribution, including
zero iterations. This gives algebraic component closure and a rate covariance
estimate from iteration-to-iteration fluctuations (only ten chunks per part).
"""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path

PAIRS={'ff_shift':('ff_matched','ff_born'),'ff_beyond_order1':('ff_matched','ff_expanded'),
       'hard_shift':('hard_projected','hard_geometry'),
       'decay_order1':('ff_expanded','hard_shift'),'decay_matched':('ff_matched','hard_shift'),
       'photos_shift':('photos_post','photos_born'),'photos_minus_ff':('photos_post','ff_matched'),
       'photos_minus_decay':('photos_post','decay_matched'),'photos_minus_order1':('photos_post','decay_order1'),
       'dressed_matched':('decay_matched','dressing_shift'),
       'photos_minus_dressed':('photos_post','dressed_matched'),
       'dressed_minus_order1':('dressed_matched','decay_order1'),
       'coupled_shift':('coupled_projected','coupled_geometry'),
       'coupled_matched':('ff_matched','coupled_shift'),
       'coupled_minus_order1':('coupled_matched','decay_order1'),
       'photos_minus_coupled':('photos_post','coupled_matched')}
SUMS={'hard_shift','decay_order1','decay_matched','dressed_matched','coupled_shift','coupled_matched'}

def analyse(path):
    meta=json.loads((path/'manifest.json').read_text())
    coupled=meta.get('environment',{}).get('LEPTONFF_COUPLED_HARD')=='1'
    dressed=meta.get('environment',{}).get('LEPTONFF_DRESS_HARD')=='1'
    groups=defaultdict(lambda:defaultdict(dict));maxclosure=0
    with (path/'iteration_rates.csv').open() as file:
        for row in csv.reader(file):
            stage,part,ips,iteration,ncall=map(int,row[:5]);label=row[5];value,variance=map(float,row[6:])
            if label.startswith('sp_'):continue
            if stage!=1:raise RuntimeError('Only production-stage rows are supported')
            if label in groups[(part,ips)][iteration]:raise RuntimeError('Duplicate iteration; possible restart')
            if variance<0 or not all(map(math.isfinite,[value,variance])):raise RuntimeError('Invalid raw iteration moments')
            groups[(part,ips)][iteration][label]=(value,variance,ncall)
    labels=sorted(next(iter(next(iter(groups.values())).values())))
    totals={key:[0,0] for key in labels};weighted={key:[0,0] for key in labels}
    covariance=[[0. for _ in labels] for _ in labels];counts={}
    for group,data in groups.items():
        iterations=list(data.values());n=len(iterations);counts[str(group)]=n
        if n<2 or any(sorted(it)!=labels for it in iterations):raise RuntimeError('Incomplete iteration records')
        if len({it[labels[0]][2] for it in iterations})!=1:raise RuntimeError('Equal-weight audit requires fixed calls per contribution')
        means={key:sum(it[key][0] for it in iterations)/n for key in labels}
        for it in iterations:
            for result,(left,right) in PAIRS.items():
                if result not in it:continue
                if ('dressed' in result) and not dressed:continue
                if ('coupled' in result) and not coupled:continue
                if left not in it or right not in it:
                    raise RuntimeError('Incomplete raw rate columns; rerun with the updated iteration-audit filter')
                residual=it[result][0]-it[left][0]-(1 if result in SUMS else -1)*it[right][0]
                maxclosure=max(maxclosure,abs(residual)*1e-3)
        if 'dressing_shift' in labels:
            for it in iterations:
                residual=it['dressing_shift'][0]-sum(it[k][0] for k in ['boundary_dressing','projected_dressing','geometry_dressing'])
                maxclosure=max(maxclosure,abs(residual)*1e-3)
        for key in labels:
            totals[key][0]+=means[key]*1e-3
            totals[key][1]+=sum(it[key][1] for it in iterations)/(n*n)*1e-6
            nonzero=[it[key] for it in iterations if it[key][0]!=0]
            if nonzero:
                inv=sum(1/r[1] for r in nonzero)
                weighted[key][0]+=sum(r[0]/r[1] for r in nonzero)/inv*1e-3
                weighted[key][1]+=1/inv*1e-6
        for i,left in enumerate(labels):
            for j,right in enumerate(labels):
                covariance[i][j]+=sum((it[left][0]-means[left])*(it[right][0]-means[right]) for it in iterations)/(n*(n-1))*1e-6
    meta=json.loads((path/'manifest.json').read_text())
    if meta['status']!='completed':raise RuntimeError('Host run incomplete')
    maximum_reconstruction=0
    for key in labels:
        totals[key][1]=math.sqrt(totals[key][1]);weighted[key][1]=math.sqrt(weighted[key][1])
        actual=meta['histograms_fb'][key][0]
        maximum_reconstruction=max(maximum_reconstruction,abs(weighted[key][0]-actual[2]*1e-3),abs(weighted[key][1]-actual[3]*1e-3))
    if maxclosure>1e-7 or maximum_reconstruction>1e-3:raise RuntimeError('Iteration closure/MCFM reconstruction failure')
    return {'equal_iteration_rates_pb':totals,'reconstructed_mcfm_rates_pb':weighted,'rate_covariance_labels':labels,
            'rate_covariance_from_iteration_fluctuations_pb2':covariance,'iterations_per_contribution':counts,
            'maximum_iteration_algebra_residual_pb':maxclosure,'maximum_mcfm_reconstruction_residual_pb':maximum_reconstruction,
            'raw_sha256':hashlib.sha256((path/'iteration_rates.csv').read_bytes()).hexdigest(),
            'scope':'Equal-weight production iterations; fixed calls per part; ten-chunk rate covariance is a diagnostic, not full cross-bin covariance'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('runs',nargs='+',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    results={str(path):analyse(path) for path in a.runs}
    report={'runs':results}
    if len(a.runs)==2:
        left,right=[results[str(p)]['equal_iteration_rates_pb'] for p in a.runs]
        report['separate_run_controls']={}
        for key in left:
            x,ex=left[key];y,ey=right[key];err=math.hypot(ex,ey)
            report['separate_run_controls'][key]={'difference_pb':x-y,'combined_error_pb':err,'pull':(x-y)/err if err else 0}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='runs'},indent=2))
if __name__=='__main__':main()
