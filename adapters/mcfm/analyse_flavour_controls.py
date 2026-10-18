#!/usr/bin/env python3
"""Audit flavour/threshold W/Z controls with grouped-event spectrum covariance."""
import argparse,csv,hashlib,json,math
from pathlib import Path
from analyse_spectra import analyse
from run_flavour_controls import OPTIONS

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def raw(path):
    with (path/'iteration_rates.csv').open() as f:
        return {(tuple(r[:5]),r[5]):float(r[6]) for r in csv.reader(f)}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    p.add_argument('--legacy-reference',type=Path,default=Path('results/spectrum_controls_v020'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checks=[];flags=[];summary={};provenance={}
    def check(label,passed,**details):checks.append({'label':label,'passed':bool(passed),**details})
    for process in ['wplus','wminus','z']:
        paths={case:a.results/f'nlo_{process}_flavours_{case}' for case in OPTIONS}
        reports={};metas={};raws={}
        for case,path in paths.items():
            meta=json.loads((path/'manifest.json').read_text());metas[case]=meta
            if meta['environment']['LEPTONFF_COUPLED_HARD']!='1' or meta['environment']['LEPTONFF_FRACTION_CELLS']!='1':
                raise RuntimeError('Wrong matching/measurement model')
            expected={name:digest(path/name) for name in ['manifest.json','iteration_rates.csv','event_covariance.csv']}
            output=a.results/f'{process}_{case}.json'
            report=json.loads(output.read_text()) if output.exists() else None
            if not report or report['provenance']!=expected or report.get('analysis_source_sha256')!=digest(Path(__file__).with_name('analyse_spectra.py')):
                report=analyse(path);output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
            reports[case]=report;raws[case]=raw(path)
            check(f'{process} {case}: grouped rates, coverage, expansion and covariance',True,
                  bin_rate_residual_pb=report['maximum_bin_rate_or_overlap_residual_pb'],
                  covariance_closure_residual_pb2=report['maximum_covariance_bin_rate_residual_pb2'])
            if meta['integration_warnings']:flags.append({'process':process,'case':case,'warnings':meta['integration_warnings']})
        if len({m['executable_sha256'] for m in metas.values()})!=1:raise RuntimeError('Mixed campaign executables')
        nominal=raws['nominal'];zero=raws['zero'];error=0.
        for (key,label),x in zero.items():
            if not label.startswith('sp_'):continue
            stem,c=label.rsplit('_',1)
            if c in ['nlo','match']:error=max(error,abs(x-zero[(key,stem+'_born')])*1e-3)
            if c in ['order2','higher']:error=max(error,abs(x)*1e-3)
            if c=='photos':error=max(error,abs(x-zero[(key,stem+'_psborn')])*1e-3)
        check(process+' zero-coupling identities',error<1e-9,maximum_residual_pb=error)
        # Identical numerical alpha(mu0) isolates flavour effects beyond O(alpha).
        for left,right,components,tolerance in [
            ('nominal','threshold_half',['born','nlo','psborn','photos'],0),
            ('nominal','threshold_double',['born','nlo','psborn','photos'],0),
            ('msbar','legacy',['born','nlo','psborn','photos'],0),
            ('nominal','grid8192',['born','nlo','psborn','photos'],0),
            ('nominal','scale_half',['born','nlo','psborn','photos'],1e-5),
            ('nominal','scale_double',['born','nlo','psborn','photos'],1e-5),
            ('nominal','slivers',['born','psborn','photos'],0)]:
            residual=max(abs(x-raws[right][key])*1e-3 for key,x in raws[left].items()
                         if key[1].startswith('sp_') and key[1].rsplit('_',1)[1] in components)
            check(f'{process} {left}/{right}: continuum/Born/PHOTOS invariance',residual<=tolerance,
                  maximum_residual_pb=residual,tolerance_pb=tolerance)
        for name in ['histograms_fb','spectra_histograms_fb']:
            check(f'{process} preset/file {name}: bitwise',metas['nominal'][name]==metas['file'][name])
        for name in ['iteration_rates.csv','event_covariance.csv']:
            check(f'{process} preset/file {name}: bitwise',digest(paths['nominal']/name)==digest(paths['file']/name))
        gold=a.legacy_reference/f'nlo_{process}_spectra_nominal'
        old=json.loads((gold/'manifest.json').read_text())
        for name in ['histograms_fb','spectra_histograms_fb']:
            check(f'{process} v0.20 legacy {name}: bitwise',old[name]==metas['legacy'][name])
        for name in ['iteration_rates.csv','event_covariance.csv']:
            check(f'{process} v0.20 legacy {name}: bitwise',digest(gold/name)==digest(paths['legacy']/name))
        alpha=metas['nominal']['flavour_matching']['alpha_input'];mu=.1056583755;me=.00051099895069
        expected=alpha/(1-alpha*2*math.log(mu/me)/(3*math.pi))
        coupling=metas['nominal']['resolved_flavour_coupling']
        check(process+' on-shell input conversion',abs(coupling['alpha_msbar_mu0']-expected)<1e-16,
              expected_alpha_msbar_mu0=expected,resolved=coupling)
        bres=max(abs(metas[c]['resolved_flavour_coupling']['beta0_integral_mu_frag']-coupling['beta0_integral_mu_frag'])
                 for c in ['threshold_half','threshold_double','zero'])
        check(process+' running logarithm independent of threshold location and alpha',bres<1e-12,maximum_residual=bres)
        base=reports['nominal'];i=base['labels'].index('coupled_matched')
        totals={case:{'matched_pb':r['mean_pb'][i],
                      'statistical_error_pb':math.sqrt(max(0.,r['covariance_pb2'][i][i])),
                      'nlo_pb':r['mean_pb'][r['labels'].index('decay_order1')],
                      'beyond_order1_pb':r['mean_pb'][r['labels'].index('coupled_minus_order1')]}
                for case,r in reports.items()}
        sensitivities={}
        for obs,rows in base['spectra'].items():
            sensitivities[obs]=[{'lower':row['lower'],'upper':row['upper'],
                 'nominal_pb':row['integrals_pb']['match'],
                 'case_shifts_pb':{case:reports[case]['spectra'][obs][k]['integrals_pb']['match'][0]-row['integrals_pb']['match'][0]
                                   for case in OPTIONS if case!='nominal'}} for k,row in enumerate(rows)]
        summary[process]={'total_rates':totals,'spectrum_sensitivities':sensitivities}
        provenance[process]={case:{'manifest_sha256':digest(path/'manifest.json'),
                                  'raw_rates_sha256':digest(path/'iteration_rates.csv'),
                                  'event_covariance_sha256':digest(path/'event_covariance.csv'),
                                  'executable_sha256':metas[case]['executable_sha256'],
                                  'source_sha256':metas[case]['source_sha256']} for case,path in paths.items()}
    failed=[c for c in checks if not c['passed']]
    report={'required_checks':len(checks),'failed_checks':len(failed),'checks':checks,
            'diagnostic_flags':flags,'summary':summary,'provenance':provenance,
            'analysis_source_sha256':digest(Path(__file__)),
            'scope':'Leptonic variable-flavour NLL final-state QED on NLO QCD W±/Z production; coloured entries are pure-QED test controls',
            'limitations':['Reduced-call validation, not differential precision certification',
                'Threshold variations probe uncomputed terms and numerical sensitivity; not a complete uncertainty budget',
                'No hadronic vacuum polarisation or QCD evolution of coloured fragmentation parents',
                'Extra-pair selection is the inclusive tagged/distinct-parent definition',
                'PHOTOS differs in running-coupling and pair content; comparison is diagnostic',
                'PDF/alpha_s, low-recoil resummation and data-definition alignment remain open']}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'required_checks':len(checks),'failed_checks':len(failed),'host_warning_flags':len(flags)},indent=2))
    if failed:raise SystemExit('Required flavour controls failed: '+str([c['label'] for c in failed]))

if __name__=='__main__':main()
