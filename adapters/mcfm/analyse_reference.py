#!/usr/bin/env python3
"""Collect the NLO W/Z interface pilot and validate identity/control runs."""
import argparse
import json
import math
from pathlib import Path

def read(directory):
    p=directory/'manifest.json';m=json.loads(p.read_text())
    if m['status']!='completed' or m['exit_code']!=0:
        raise RuntimeError('Incomplete/failed run: '+str(directory))
    return m

def rates(m):
    return {name:{'value_pb':rows[0][2]*.001,'error_pb':rows[0][3]*.001}
      for name,rows in m['histograms_fb'].items()}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=Path('results'))
    parser.add_argument('--output',type=Path,default=Path('results/nlo_host_reference.json'))
    args=parser.parse_args();root=args.results
    samples={}
    for process in ['wplus','wminus','z']:
        samples[process]=read(root/('nlo_'+process+'_photos_seeded'))
    summary={'scope':'NLO QCD plus one-family Born-parent FF-only pilot; finite W/Z hard O(alpha) absent',
      'photos':'3.56; exponentiation/interference on; W ME on, Z ME off; alpha(0)=1/137.035999084',
      'rates':{name:rates(m) for name,m in samples.items()},'checks':{},
      'warning':'Printed histogram means use independent per-bin iteration weighting. Use directly accumulated difference histograms for paired errors; no exact closure of printed means or full cross-bin covariance is asserted.'}
    for process in ['wplus','wminus','z']:
        zero=read(root/('nlo_'+process+'_zero'));zr=rates(zero)
        b=zr['ff_born']['value_pb']
        closure=max(abs(zr[k]['value_pb']-b) for k in ['ff_expanded','ff_resummed','ff_matched'])
        if closure>1e-9 or abs(zr['ff_shift']['value_pb'])>1e-9:
            raise RuntimeError('Zero-alpha identity failed')
        summary['checks'][process+'_zero_alpha_max_residual_pb']=closure
    for process in ['wplus','wminus','z']:
        nr=rates(read(root/('nlo_'+process+'_native')))['total_cross']
        br=summary['rates'][process]['ff_born']
        pull=(br['value_pb']-nr['value_pb'])/math.hypot(br['error_pb'],nr['error_pb'])
        if abs(pull)>4:raise RuntimeError('Native identity comparison outside four combined standard errors')
        summary['checks'][process+'_native_identity_pull']=pull
    for process,fine_name in [('wplus','nlo_wplus_grid8192'),('wminus','nlo_wminus_grid8192'),('z','nlo_z_grid4096')]:
        base=summary['rates'][process];fine=rates(read(root/fine_name))
        # A changed FF grid must not change the underlying QCD or PHOTOS stream.
        for key in ['total_cross','ff_born','photos_post']:
            if base[key]!=fine[key]:raise RuntimeError('Grid scan did not retain common production/PHOTOS')
        summary['checks'][process+'_grid_matched_shift_pb']=fine['ff_matched']['value_pb']-base['ff_matched']['value_pb']
        varied=rates(read(root/('nlo_'+process+'_alpha01')))
        for key in ['ff_born','ff_matched','photos_post','photos_minus_ff']:
            pull=(varied[key]['value_pb']-base[key]['value_pb'])/math.hypot(varied[key]['error_pb'],base[key]['error_pb'])
            if abs(pull)>4:raise RuntimeError('Dipole-alpha stability outside four combined standard errors')
            summary['checks'][process+'_dipole_alpha_'+key+'_pull']=pull
    summary['run_manifests']={str(p.relative_to(root)):json.loads(p.read_text())
      for p in root.glob('nlo_*/manifest.json') if p.parent.name in {
        'nlo_'+process+suffix for process in ['wplus','wminus','z'] for suffix in
        ['_photos_seeded','_native','_zero','_alpha01','_grid8192','_grid4096']}}
    args.output.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({'checks':summary['checks'],'output':str(args.output)},indent=2))

if __name__=='__main__':main()
