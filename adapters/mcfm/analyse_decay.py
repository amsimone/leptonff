#!/usr/bin/env python3
"""Audit the optional finite-decay MCFM controls; no fitted hard corrections."""
import argparse,hashlib,json,math
from pathlib import Path

PROCESSES=['wplus','wminus','z']

def load(root,process,tag):
    path=root/f'nlo_{process}_decay_{tag}'/'manifest.json'
    x=json.loads(path.read_text())
    if x['status']!='completed' or x['exit_code']!=0:raise RuntimeError(f'Incomplete run: {path}')
    if x.get('integration_warnings'):raise RuntimeError(f'Integration warning: {path}')
    if hashlib.sha256((path.parent/'input.ini').read_bytes()).hexdigest()!=x['config_sha256']:
        raise RuntimeError(f'Configuration hash mismatch: {path}')
    return x

def rate(x,key):
    r=x['histograms_fb'][key]
    if len(r)!=1 or len(r[0])!=4 or not all(math.isfinite(v) for v in r[0]):raise RuntimeError('Invalid rate histogram')
    return [1e-3*r[0][2],1e-3*r[0][3]]

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,default=Path('results'))
    p.add_argument('--output',type=Path,default=Path('results/decay_host_reference.json'));a=p.parse_args()
    checks=[];rows={};provenance={}
    def exact(left,right,keys,label):
        if any(rate(left,k)!=rate(right,k) for k in keys):raise RuntimeError('Exact control failed: '+label)
        checks.append({'label':label,'kind':'exact','passed':True})
    def compare(left,right,keys,label,required=True):
        values={}
        for k in keys:
            x,ex=rate(left,k);y,ey=rate(right,k);den=math.hypot(ex,ey)
            pull=(x-y)/den if den else 0
            if not den and x!=y:raise RuntimeError('Zero-error shift: '+label)
            values[k]={'shift_pb':x-y,'combined_error_pb':den,'pull':pull}
            if required and abs(pull)>3:raise RuntimeError(f'Control exceeds 3 combined sigma: {label} {k} {pull}')
        checks.append({'label':label,'kind':'separate-run combined errors (cross-run covariance not stored)',
                       'values':values,'passed':all(abs(v['pull'])<=3 for v in values.values()),'required':required})
    for proc in PROCESSES:
        runs={tag:load(a.results,proc,tag) for tag in ['order1','exponentiated','zero','slivers','grid8192','alpha01']}
        x=runs['exponentiated'];fo=runs['order1'];zero=runs['zero'];grid=runs['grid8192'];fine=runs['slivers']
        provenance[proc]={tag:{k:r[k] for k in ['source_sha256','executable_sha256','config_sha256','environment']}
                          for tag,r in runs.items()}
        for key in ['hard_projected','hard_geometry','hard_shift']:
            if rate(zero,key)!=[0,0]:raise RuntimeError('Nonzero hard term at alpha=0')
        for key in ['ff_expanded','ff_resummed','ff_matched','decay_order1','decay_matched']:
            if rate(zero,key)!=rate(zero,'ff_born'):raise RuntimeError('Zero-alpha identity failure')
        if rate(zero,'photos_post')!=rate(zero,'photos_born') or rate(zero,'photos_shift')!=[0,0]:
            raise RuntimeError('Zero-alpha PHOTOS identity failure')
        checks.append({'label':proc+' zero-coupling hard/FF/PHOTOS identities','kind':'exact','passed':True})
        exact(x,fo,['ff_born','ff_expanded','ff_resummed','ff_matched','hard_projected','hard_geometry','hard_shift',
                    'decay_order1','decay_matched'],proc+' analytic branch independent of generator mode')
        exact(x,grid,['ff_born','ff_expanded','hard_projected','hard_geometry','hard_shift','decay_order1',
                      'photos_born','photos_post'],proc+' grid changes only resummed FF branch')
        exact(x,fine,['ff_born','ff_expanded','ff_matched','hard_projected','photos_born','photos_post'],
              proc+' geometric exclusions change only geometric branch')
        compare(fine,x,['hard_geometry','hard_shift','decay_matched','photos_minus_decay'],proc+' geometric cutoff refinement')
        compare(runs['alpha01'],x,['ff_born','hard_shift','decay_matched',
                                 'photos_minus_decay'],proc+' QCD dipole alpha 1 to 0.1 physical sums')
        compare(runs['alpha01'],x,['hard_projected','hard_geometry'],
                proc+' QCD dipole component diagnostic',required=False)
        if proc=='wminus':
            refined=[load(a.results,proc,tag) for tag in ['precision1','precision01']]
            compare(refined[0],refined[1],['ff_born','hard_projected','hard_geometry','hard_shift',
                                        'decay_matched','photos_minus_decay'],
                    'W- independent seed 2089, four-times-call subtraction control')
            for tag,r in zip(['precision1','precision01'],refined):
                provenance[proc][tag]={k:r[k] for k in ['source_sha256','executable_sha256','config_sha256','environment']}
        rows[proc]={k:rate(x,k) for k in ['ff_born','ff_matched','hard_projected','hard_geometry','hard_shift',
                      'decay_order1','decay_matched','photos_post','photos_minus_decay']}
        rows[proc]['photos_order1_minus_analytic_order1']=rate(fo,'photos_minus_order1')
        rows[proc]['grid_refinement_shift_pb']=rate(grid,'decay_matched')[0]-rate(x,'decay_matched')[0]
        rows[proc]['maximum_chisq_per_iteration']=x['maximum_chisq_per_iteration']
    iteration_path=a.results/'wminus_iteration_comparison.json'
    iteration=json.loads(iteration_path.read_text())
    for tag,base in [('audit1','precision1'),('audit01','precision01')]:
        instrumented=load(a.results,'wminus',tag);original=load(a.results,'wminus',base)
        if instrumented['histograms_fb']!=original['histograms_fb']:
            raise RuntimeError('Instrumentation changed host rate output')
        checks.append({'label':'W- '+tag+' instrumentation leaves histograms bitwise identical',
                       'kind':'exact','passed':True})
    controls=iteration['separate_run_controls']
    selected=['hard_projected','hard_geometry','hard_shift','decay_matched','photos_minus_decay']
    if any(abs(controls[key]['pull'])>3 for key in selected):
        raise RuntimeError('Common-iteration W- subtraction control failed')
    checks.append({'label':'W- common-iteration subtraction control',
                   'kind':'equal iteration weights, independent of each histogram value/variance',
                   'values':{key:controls[key] for key in selected},'passed':True})
    report={'scope':'NLO QCD production times unitary photonic additive NLO decay + LL; integration errors only, not hard-dressed NLL',
            'iteration_audit':str(iteration_path),'validation_status':'physical-rate controls pass; original projected-component statistical diagnostic retained',
            'rates_pb':rows,'checks':checks,'component_flags':[c for c in checks if not c['passed']],'provenance':provenance,
            'statistics_note':'MCFM combines iterations independently per histogram. Paired differences are accumulated directly; do not subtract independently averaged rate rows.'}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'controls_passed':sum(c['passed'] for c in checks),'diagnostic_flags':len(report['component_flags']),'rates_pb':rows},indent=2))
if __name__=='__main__':main()
