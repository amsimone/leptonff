#!/usr/bin/env python3
"""Audit the complete spectrum campaign and separate uncertainty sources."""
import argparse,csv,hashlib,json,math
from pathlib import Path
from analyse_spectra import analyse,COMPONENTS,definitions
from run_spectrum_controls import OPTIONS

def raw(path):
    with (path/'iteration_rates.csv').open() as file:
        return {(tuple(row[:5]),row[5]):float(row[6]) for row in csv.reader(file) if row[5].startswith('sp_')}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checks=[];flags=[];summary={};provenance={}
    def check(label,passed,**details):checks.append({'label':label,'passed':passed,**details})
    for process in ['wplus','wminus','z']:
        paths={case:a.results/f'nlo_{process}_spectra_{case}' for case in OPTIONS}
        reports={};meta={};raws={}
        for case,path in paths.items():
            m=json.loads((path/'manifest.json').read_text());meta[case]=m
            if hashlib.sha256((path/'input.ini').read_bytes()).hexdigest()!=m['config_sha256']:raise RuntimeError('Changed configuration')
            if m['environment']['LEPTONFF_FRACTION_CELLS']!='1' or m['environment']['LEPTONFF_COUPLED_HARD']!='1':raise RuntimeError('Wrong theory branch')
            output=a.results/f'{process}_{case}.json'
            expected={name:hashlib.sha256((path/name).read_bytes()).hexdigest() for name in ['manifest.json','iteration_rates.csv','event_covariance.csv']}
            r=json.loads(output.read_text()) if output.exists() else None
            if not r or r['provenance']!=expected or r.get('analysis_source_sha256')!=hashlib.sha256(Path(__file__).with_name('analyse_spectra.py').read_bytes()).hexdigest():
                r=analyse(path);output.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
            reports[case]=r;raws[case]=raw(path)
            check(f'{process} {case}: grouped sums, covariance, coverage and overlap',True,
                  maximum_bin_residual_pb=r['maximum_bin_rate_or_overlap_residual_pb'],
                  maximum_covariance_residual_pb2=r['maximum_covariance_bin_rate_residual_pb2'])
            if m['integration_warnings']:flags.append({'process':process,'case':case,'host_warnings':m['integration_warnings']})
        if len({m['executable_sha256'] for m in meta.values()})!=1:raise RuntimeError('Campaign mixes adapter executables')
        base=reports['nominal'];zero=reports['zero'];zraw=raws['zero'];zerr=0.
        for (iteration,label),x in zraw.items():
            stem,component=label.rsplit('_',1)
            if component in ['nlo','match']:zerr=max(zerr,abs(x-zraw[(iteration,stem+'_born')])*1e-3)
            if component in ['higher','order2']:zerr=max(zerr,abs(x)*1e-3)
            if component=='photos':zerr=max(zerr,abs(x-zraw[(iteration,stem+'_psborn')])*1e-3)
            if component=='psdiff':zerr=max(zerr,abs(x-zraw[(iteration,stem+'_psborn')]+zraw[(iteration,stem+'_born')])*1e-3)
        check(process+' zero coupling with separate massive PHOTOS Born map',zerr<1e-9,maximum_residual_pb=zerr)
        for case,components,tolerance in [('grid8192',['born','nlo','psborn','photos'],0),
                                          ('scale_half',['born','nlo','psborn','photos'],1e-5),
                                          ('scale_double',['born','nlo','psborn','photos'],1e-5),
                                          ('slivers',['born','psborn','photos'],0)]:
            residual=max(abs(x-raws[case][key])*1e-3 for key,x in raws['nominal'].items() if key[1].rsplit('_',1)[1] in components)
            check(process+' '+case+' continuum/Born/PHOTOS invariance',residual<=tolerance,maximum_residual_pb=residual,tolerance_pb=tolerance)
        budget={}
        for observable,rows in base['spectra'].items():
            bins=[]
            for k,row in enumerate(rows):
                central,error=row['integrals_pb']['match']
                entry={'lower':row['lower'],'upper':row['upper'],'central_pb':central,'statistical_error_pb':error}
                for name,cases in [('qcd_seven_point',['nominal','qcd_half','qcd_double','qcd_rhalf','qcd_rdouble','qcd_fhalf','qcd_fdouble']),
                                   ('fragmentation_three_point',['nominal','scale_half','scale_double'])]:
                    points={c:reports[c]['spectra'][observable][k]['integrals_pb']['match'] for c in cases}
                    entry[name]={'points_pb':points,'envelope_shift_pb':[min(v[0] for v in points.values())-central,max(v[0] for v in points.values())-central]}
                for case in ['grid8192','slivers','seed2089']:
                    value,err=reports[case]['spectra'][observable][k]['integrals_pb']['match']
                    diagnostic_error=math.hypot(err,error)
                    entry[case]={'shift_pb':value-central,'unpaired_quadrature_error_pb':diagnostic_error,
                                 'diagnostic_pull':(value-central)/diagnostic_error if diagnostic_error else 0}
                entry['photos_minus_match_pb']=row['integrals_pb']['psdiff']
                bins.append(entry)
            budget[observable]=bins
        summary[process]={'nominal_total_rates_pb':{label:[base['mean_pb'][i],math.sqrt(max(0.,base['covariance_pb2'][i][i]))] for i,label in enumerate(base['labels'][:8])},
                          'per_bin_uncertainty_components':budget}
        provenance[process]={case:{k:m[k] for k in ['source_sha256','executable_sha256','config_sha256','requested_inputs','environment']} for case,m in meta.items()}
    failures=[c for c in checks if not c['passed']]
    report={'checks':checks,'required_checks':len(checks),'failed_checks':len(failures),'diagnostic_flags':flags,
            'summary':summary,'provenance':provenance,
            'scope':'Physical observed W and charge-symmetric Z spectra; full grouped-event covariance and scale/refinement diagnostics in the one-family MSbar model',
            'limitations':['Reduced-call controls test functionality and sensitivity; scale envelopes contain Monte Carlo noise',
                           'Grid and cutoff shifts are sensitivities, not certified numerical errors',
                           'No cross-run event products: refinement pulls use unpaired quadrature diagnostics',
                           'Physical flavours, threshold/input-scheme matching and PDF uncertainty remain open',
                           'PHOTOS has different pair/coupling content; the coupled comparison is diagnostic',
                           'Fixed-order NLO QCD recoil is unreliable for precision Jacobian/low-recoil shapes without resummation']}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'required_checks':len(checks),'failed_checks':len(failures),'host_warning_flags':len(flags)},indent=2))
    if failures:raise SystemExit('Required spectrum controls failed')
if __name__=='__main__':main()
