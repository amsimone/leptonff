#!/usr/bin/env python3
"""Audit HVP/LO-QCD pointlike benchmark, keeping approximation and host flags."""
import argparse,csv,hashlib,json,math,sys
from pathlib import Path
from analyse_spectra import analyse
from run_hadronic_controls import OPTIONS
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'python'))
from check_flavour_reference import running

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def raw(path):
    with (path/'iteration_rates.csv').open() as f:return {(tuple(r[:5]),r[5]):float(r[6]) for r in csv.reader(f)}
def expected_anchor(meta,path):
    data={}
    for line in (path/'hadronic.dat').read_text().splitlines():
        words=line.split('#',1)[0].split()
        if words:data[words[0]]=float(words[1])
    alpha=meta['flavour_matching']['alpha_input'];pi=math.pi
    a=alpha/(1-alpha*2*math.log(.1056583755/.00051099895069)/(3*pi))/(2*pi)
    a=running(a,2*math.log(1.77693/.1056583755),4/3,2)
    a=running(a,2*math.log(data['reference_scale']/1.77693),2,3)
    delta=data['delta_had5_os']+meta['hadronic_input']['input_standard_deviation_shift']*data['uncertainty']
    strong=meta['hadronic_input']['alpha_s_override']
    if strong is None:strong=data['alpha_s_reference']
    shift=delta*alpha/(1/137.035999084)+alpha/pi*(55/27+strong/pi*(605/108-44*1.2020569031595942854/9))
    return 1/(1/(2*pi*a)-shift/alpha) if alpha else 0.

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    p.add_argument('--legacy-reference',type=Path,default=Path('results/flavour_matching_v021'))
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    checks=[];flags=[];summary={};provenance={}
    def check(label,passed,**details):checks.append({'label':label,'passed':bool(passed),**details})
    for process in ['wplus','wminus','z']:
        paths={case:a.results/f'nlo_{process}_hadronic_{case}' for case in OPTIONS};metas={};reports={};raws={}
        for case,path in paths.items():
            m=json.loads((path/'manifest.json').read_text());metas[case]=m
            out=a.results/f'{process}_{case}.json';r=json.loads(out.read_text()) if out.exists() else None
            provenance_now={name:digest(path/name) for name in ['manifest.json','iteration_rates.csv','event_covariance.csv']}
            if not r or r['provenance']!=provenance_now or r['analysis_source_sha256']!=digest(Path(__file__).with_name('analyse_spectra.py')):
                r=analyse(path);out.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
            reports[case]=r;raws[case]=raw(path)
            check(f'{process} {case}: rates, spectra, overlap, grouped covariance',True,
                  bin_rate_residual_pb=r['maximum_bin_rate_or_overlap_residual_pb'],covariance_closure_pb2=r['maximum_covariance_bin_rate_residual_pb2'])
            if m['integration_warnings']:flags.append({'process':process,'case':case,'warnings':m['integration_warnings']})
        if len({m['executable_sha256'] for m in metas.values()})!=1:raise RuntimeError('Mixed executables')
        for case in OPTIONS:
            if case in ['nominal','zero','legacy']:continue
            residual=max(abs(x-raws[case][k])*1e-3 for k,x in raws['nominal'].items()
                         if k[1].startswith('sp_') and k[1].rsplit('_',1)[1] in ['born','nlo','psborn','photos'])
            check(f'{process} {case}: fixed-order continuum/Born/PHOTOS invariant',residual<=1e-5,maximum_residual_pb=residual)
        zero=raws['zero'];error=0.
        for (it,label),value in zero.items():
            if not label.startswith('sp_'):continue
            stem,c=label.rsplit('_',1)
            if c in ['nlo','match']:error=max(error,abs(value-zero[(it,stem+'_born')])*1e-3)
            if c in ['order2','higher']:error=max(error,abs(value)*1e-3)
            if c=='photos':error=max(error,abs(value-zero[(it,stem+'_psborn')])*1e-3)
        check(process+' zero alpha identities with nonzero strong coupling',error<1e-9,maximum_residual_pb=error)
        old=a.legacy_reference/f'nlo_{process}_flavours_legacy';gold=json.loads((old/'manifest.json').read_text())
        for name in ['histograms_fb','spectra_histograms_fb']:
            check(process+' legacy '+name+' bitwise',metas['legacy'][name]==gold[name])
        for name in ['iteration_rates.csv','event_covariance.csv']:
            check(process+' legacy '+name+' bitwise',digest(paths['legacy']/name)==digest(old/name))
        central=metas['nominal']['resolved_flavour_coupling']
        check(process+' turning QCD off preserves the HVP coupling input',central==metas['qcd_off']['resolved_flavour_coupling'])
        residual=max(abs(metas[c]['resolved_flavour_coupling']['alpha_msbar_mu_frag']-central['alpha_msbar_mu_frag']) for c in ['match_low','match_high'])
        check(process+' matching-scale changes preserve the high-scale anchored coupling',residual<2e-14,maximum_residual=residual)
        check(process+' input uncertainty moves alpha monotonically',
              metas['had_minus']['resolved_flavour_coupling']['alpha_msbar_mu_frag']<central['alpha_msbar_mu_frag']<metas['had_plus']['resolved_flavour_coupling']['alpha_msbar_mu_frag'])
        if process=='z':
            for case in ['nominal','qcd_off','had_minus','had_plus','match_low','match_high','grid1024','alternative']:
                residual=abs(metas[case]['resolved_flavour_coupling']['alpha_msbar_mu_frag']-expected_anchor(metas[case],paths[case]))
                check('Z '+case+' sourced high-scale anchor',residual<1e-14,maximum_residual=residual)
        base=reports['nominal'];totals={}
        for case,r in reports.items():
            totals[case]={label:[r['mean_pb'][i],math.sqrt(max(0.,r['covariance_pb2'][i][i]))] for i,label in enumerate(r['labels'][:8])}
        shifts={obs:[{'lower':row['lower'],'upper':row['upper'],'nominal_pb':row['integrals_pb']['match'],
            'shifts_pb':{case:reports[case]['spectra'][obs][k]['integrals_pb']['match'][0]-row['integrals_pb']['match'][0] for case in OPTIONS if case!='nominal'}}
            for k,row in enumerate(rows)] for obs,rows in base['spectra'].items()}
        summary[process]={'total_rates_pb':totals,'spectrum_sensitivities':shifts}
        provenance[process]={case:{'manifest_sha256':digest(path/'manifest.json'),'raw_rates_sha256':digest(path/'iteration_rates.csv'),
            'event_covariance_sha256':digest(path/'event_covariance.csv'),'executable_sha256':metas[case]['executable_sha256'],'source_sha256':metas[case]['source_sha256']} for case,path in paths.items()}
    failed=[c for c in checks if not c['passed']]
    report={'required_checks':len(checks),'failed_checks':len(failed),'checks':checks,'diagnostic_flags':flags,
        'summary':summary,'provenance':provenance,'analysis_source_sha256':digest(Path(__file__)),
        'scope':'Sourced hadronic alpha input and LO QCD coloured-parent evolution with an explicit pointlike fragmentation boundary on NLO QCD W/Z production',
        'limitations':['Not a fitted nonperturbative hadronic fragmentation boundary','Mixed alpha*alpha_s splitting kernels and beta terms absent',
            'Hadronic OS/MSbar conversion truncated through alpha*alpha_s; mass power terms and W/top loops omitted',
            'Input uncertainty is one contribution, not a complete phenomenology error',
            'Reduced-call and grid/matching variations are functional sensitivities; no cross-run covariance or significance claims',
            'Production uses the existing NNLO PDF reference; PDF/alpha_s and recoil/data-alignment work remain',
            'PHOTOS pair content and fixed coupling differ']}
    a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'required_checks':len(checks),'failed_checks':len(failed),'host_flags':len(flags)},indent=2))
    if failed:raise SystemExit('Hadronic controls failed: '+str([c['label'] for c in failed]))
if __name__=='__main__':main()
