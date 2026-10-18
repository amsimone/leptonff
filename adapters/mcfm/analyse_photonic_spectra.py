#!/usr/bin/env python3
"""Audit the coupling/pair-aligned PHOTOS versus photonic spectrum pilot."""
import argparse,hashlib,json
from pathlib import Path
from analyse_spectra import analyse
from analyse_spectrum_controls import raw

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    p.add_argument('--coupled-results',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checks=[];flags=[];summary={}
    def check(label,residual,tolerance=1e-9):
        checks.append({'label':label,'maximum_residual_pb':residual,'tolerance_pb':tolerance,'passed':residual<=tolerance})
    for process in ['wplus','wminus','z']:
        reports={};raws={};metas={}
        for case in ['nominal','zero','grid8192','slivers']:
            path=a.results/f'nlo_{process}_spectra_{case}';meta=json.loads((path/'manifest.json').read_text());metas[case]=meta
            if meta['environment']['LEPTONFF_PHOTONIC_LL']!='1' or meta['environment']['LEPTONFF_DRESS_HARD']!='1':raise RuntimeError('Wrong photonic branch')
            if meta['environment']['LEPTONFF_COUPLED_HARD']!='0':raise RuntimeError('Coupled branch present')
            if hashlib.sha256((path/'input.ini').read_bytes()).hexdigest()!=meta['config_sha256']:raise RuntimeError('Configuration changed')
            r=analyse(path);reports[case]=r;raws[case]=raw(path)
            (a.results/f'{process}_{case}.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
            check(process+' '+case+' grouped covariance and bin/rate closure',r['maximum_bin_rate_or_overlap_residual_pb'])
            if meta['integration_warnings']:flags.append({'process':process,'case':case,'warnings':meta['integration_warnings']})
        zero=raws['zero'];residual=0.
        for (it,label),x in zero.items():
            stem,c=label.rsplit('_',1)
            if c in ['nlo','match']:residual=max(residual,abs(x-zero[it,stem+'_born'])*1e-3)
            if c in ['ff2','higher']:residual=max(residual,abs(x)*1e-3)
            if c=='photos':residual=max(residual,abs(x-zero[it,stem+'_psborn'])*1e-3)
        check(process+' zero coupling',residual)
        for case,components in [('grid8192',['born','nlo','psborn','photos']),('slivers',['born','psborn','photos'])]:
            residual=max(abs(x-raws[case][key])*1e-3 for key,x in raws['nominal'].items() if key[1].rsplit('_',1)[1] in components)
            check(process+' '+case+' invariant components',residual,0)
        coupled=raw(a.coupled_results/f'nlo_{process}_spectra_nominal')
        keys=[key for key in raws['nominal'] if key[1].rsplit('_',1)[1] in ['born','nlo','psborn','photos']]
        if any(key not in coupled for key in keys):raise RuntimeError('Coupled reference call budget/seed differs')
        residual=max(abs(raws['nominal'][key]-coupled[key])*1e-3 for key in keys)
        check(process+' independent LL/NLL backends agree at continuum NLO and PHOTOS',residual,1e-7)
        nominal=reports['nominal']
        summary[process]={'total_rates_pb':{label:[nominal['mean_pb'][i],nominal['covariance_pb2'][i][i]**.5] for i,label in enumerate(nominal['labels'][:8])},
                          'comparison':'Fixed alpha, no pairs/running in both tools; analytic LL with at most one finite insertion; PHOTOS exponentiation, W ME on and Z ME off; leading-power vs massive kinematics',
                          'ff2_column':'FF-only LL expansion coefficient, not full dressed order2'}
    report={'required_checks':len(checks),'failed_checks':sum(not c['passed'] for c in checks),
            'checks':checks,'diagnostic_flags':flags,'summary':summary,
            'scope':'Aligned coupling/pair-content benchmark; numerical and statistical pilot, not a physical-flavour NLL prediction'}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n');print(json.dumps({k:report[k] for k in ['required_checks','failed_checks']},indent=2))
    if report['failed_checks']:raise SystemExit('Required photonic spectrum controls failed')
if __name__=='__main__':main()
