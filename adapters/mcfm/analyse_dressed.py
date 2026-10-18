#!/usr/bin/env python3
"""Audit LL dressing on the same NLO QCD contributions and raw iterations."""
import argparse,json,math,hashlib
from pathlib import Path
from analyse_iterations import analyse

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,default=Path('results'))
    p.add_argument('--output',type=Path,default=Path('results/photonic_dressing_host.json'));a=p.parse_args()
    checks=[];flags=[];rates={};allruns={};provenance={}
    def exact(left,right,keys,label):
        if any(left['histograms_fb'][k]!=right['histograms_fb'][k] for k in keys):raise RuntimeError('Exact control failed: '+label)
        checks.append({'label':label,'kind':'same-stream exact','passed':True})
    def comparison(left,right,keys,label,required=True):
        controls={}
        for key in keys:
            x,ex=left['equal_iteration_rates_pb'][key];y,ey=right['equal_iteration_rates_pb'][key]
            error=math.hypot(ex,ey);pull=(x-y)/error if error else 0
            if not error and x!=y:raise RuntimeError('Zero-error shift: '+label)
            controls[key]={'shift_pb':x-y,'combined_error_pb':error,'pull':pull}
        passed=all(abs(v['pull'])<=3 for v in controls.values())
        row={'label':label,'kind':'common iteration weights; separate-run combined errors, no cross-run covariance',
             'values':controls,'passed':passed,'required':required}
        checks.append(row)
        if not passed:
            flags.append(row)
            if required:raise RuntimeError('Required control failed: '+label)
    for process in ['wplus','wminus','z']:
        paths={tag:a.results/f'nlo_{process}_dressed_{tag}' for tag in ['audit','zero','grid8192','alpha01','slivers']}
        meta={tag:json.loads((path/'manifest.json').read_text()) for tag,path in paths.items()}
        for tag,row in meta.items():
            if row['status']!='completed' or row['exit_code']!=0 or row['integration_warnings']:
                raise RuntimeError('Incomplete/warned host run: '+str(paths[tag]))
            if hashlib.sha256((paths[tag]/'input.ini').read_bytes()).hexdigest()!=row['config_sha256']:
                raise RuntimeError('Config hash differs')
        raw={tag:analyse(path) for tag,path in paths.items()};allruns[process]=raw
        nominal=meta['audit'];zero=meta['zero'];grid=meta['grid8192'];fine=meta['slivers']
        old=json.loads((a.results/f'nlo_{process}_decay_exponentiated/manifest.json').read_text())
        exact(nominal,old,list(old['histograms_fb']),process+' preserves all v0.16 baseline columns')
        for key in ['boundary_dressing','projected_dressing','geometry_dressing','dressing_shift']:
            if zero['histograms_fb'][key][0][2:]!=[0,0]:raise RuntimeError('Nonzero dressing at zero alpha')
        checks.append({'label':process+' zero coupling dressing vanishes','kind':'exact identity','passed':True})
        for key in ['decay_order1','decay_matched','dressed_matched']:
            if zero['histograms_fb'][key]!=zero['histograms_fb']['ff_born']:raise RuntimeError('Zero alpha identity failed')
        exact(nominal,grid,['ff_born','ff_expanded','hard_projected','hard_geometry','hard_shift','decay_order1','photos_post'],
              process+' FF grid leaves continuum NLO and generator unchanged')
        exact(nominal,fine,['ff_born','ff_expanded','ff_matched','hard_projected','boundary_dressing','projected_dressing','photos_post'],
              process+' slivers affect only geometric branch')
        comparison(raw['alpha01'],raw['audit'],['dressed_matched','dressing_shift','photos_minus_dressed'],process+' QCD dipole alpha')
        comparison(raw['slivers'],raw['audit'],['geometry_dressing','dressing_shift','dressed_matched','photos_minus_dressed'],process+' finite geometry refinement')
        comparison(raw['alpha01'],raw['audit'],['boundary_dressing','projected_dressing','geometry_dressing'],process+' component diagnostic',required=False)
        rates[process]={key:raw['audit']['equal_iteration_rates_pb'][key] for key in ['boundary_dressing','projected_dressing','geometry_dressing',
                           'dressing_shift','decay_matched','dressed_matched','photos_minus_dressed','dressed_minus_order1']}
        rates[process]['grid_refinement_shift_pb']={key:raw['grid8192']['equal_iteration_rates_pb'][key][0]-raw['audit']['equal_iteration_rates_pb'][key][0]
                                                  for key in ['dressing_shift','dressed_matched','geometry_dressing']}
        provenance[process]={tag:{key:row[key] for key in ['source_sha256','executable_sha256','config_sha256','environment']}
                             for tag,row in meta.items()}
    report={'scope':'Fixed-alpha photonic LL with one dressed boundary or H1 globally; NLO QCD host; not full NLL',
            'rates_pb':rates,'checks':checks,'component_flags':flags,'raw_iteration_audits':allruns,'provenance':provenance,
            'uncertainty_scope':'Within-iteration errors combined with common iteration means; ten-chunk covariance diagnostic, no heavy-tail or cross-bin bound'}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'checks_passed':sum(c['passed'] for c in checks),'component_flags':len(flags),'rates_pb':rates},indent=2))
if __name__=='__main__':main()
