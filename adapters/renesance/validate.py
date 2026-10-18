#!/usr/bin/env python3
"""Reproduce diagnostics for the documented local Born-FSR campaign layout."""
import argparse
import json
from pathlib import Path
import numpy as np
from compare import combine, ratio


def pulls(left,right,coeff):
    values=[]
    for a,b in zip(left['bins'],right['bins']):
        x,y=ratio(a,coeff),ratio(b,coeff)
        if x and y and x['error']+y['error']:
            values.append((x['correction']-y['correction'])/np.hypot(x['error'],y['error']))
    return {'bins':len(values),'pull_rms':float(np.sqrt(np.mean(np.square(values)))),
            'max_abs_pull':max(map(abs,values))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=Path('results/born-fsr'))
    a=p.parse_args();root=a.root
    def get(folder,pattern):
        return combine([json.loads(p.read_text()) for p in (root/folder).glob(pattern)])
    fo={0:1,1:1,2:1};nll={7:1}
    fine=get('z-refinement','z-fo-*.json');ff=get('z-refinement','z-analytic-*.json')
    coarse=get('complete','z-fo-*.json');coarse_ff=get('complete','z-analytic-*.json')
    def shift(left,right,lcoeff,rcoeff):
        return max(abs(ratio(x,lcoeff)['correction']-ratio(y,rcoeff)['correction'])
                   for x,y in zip(left['bins'][2:],right['bins'][2:]))
    report={'schema':'leptonff-born-fsr-validation-v1',
        'z_order_alpha_max_difference_over_born':shift(fine,ff,fo,{5:1}),
        'z_fo_max_refinement_shift':shift(coarse,fine,fo,fo),
        'z_nll_max_refinement_shift':shift(coarse_ff,ff,nll,nll),
        'z_direct_mc_vs_quadrature':pulls(get('z-mc-control','z-fo-*.json'),fine,fo)}
    for process in ('wplus','wminus'):
        main=combine([json.loads(p.read_text()) for folder in ('central','complete')
                      for p in (root/folder).glob(process+'-fo-*.json')])
        report[process+'_soft_cutoff']=pulls(main,get('w-cutoff-control',process+'-fo-*.json'),fo)
        report[process+'_inclusive_closure']=ratio(main['bins'][0],fo)
    for process in ('z','wplus','wminus'):
        original=get('central',process+'-showers-*.json')
        if (root/'shower-cutoff-control').exists():
            control=get('shower-cutoff-control',process+'-showers-*.json')
            report[process+'_shower_cutoffs']={name:pulls(original,control,coeff)
                for name,coeff in [('photos',{3:1}),('pythia',{4:1})]}
        if (root/'pythia-mec-off').exists():
            off=get('pythia-mec-off',process+'-showers-*.json')
            report[process+'_pythia_mec']={'on':ratio(original['bins'][1],{4:1}),
                                         'off':ratio(off['bins'][1],{4:1})}
    (root/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
