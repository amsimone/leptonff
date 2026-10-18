#!/usr/bin/env python3
"""Compare conservative-cell W/Z measurements with the frozen atomic pilot."""
import argparse,hashlib,json,math,re
from pathlib import Path
from analyse_iterations import analyse
from analyse_coupled import paired_chunks

KEYS=['coupled_matched','coupled_minus_order1','coupled_order2','photos_minus_coupled']
BASE=['ff_born','ff_expanded','hard_projected','hard_geometry','hard_shift','decay_order1',
      'photos_born','photos_post','photos_shift','photos_radiated','photos_minus_order1']

def cpu(path):
    found=re.findall(r'CPU time total:\s+([\d.Ee+\-]+)',(path/'run.log').read_text())
    if not found:raise RuntimeError('Host CPU time absent')
    return float(found[-1])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=Path('results'))
    p.add_argument('--output',type=Path,default=Path('results/fraction_cells_host.json'))
    p.add_argument('--reference-results',type=Path,help='Frozen v0.18 run directory; defaults to --results')
    p.add_argument('--processes',nargs='+',choices=['wplus','wminus','z'],default=['wplus','wminus','z'])
    a=p.parse_args();reference=a.reference_results or a.results;checks=[];flags=[];raws={};summary={};provenance={}
    def exact(left,right,keys,label):
        checks.append({'label':label,'kind':'same-stream bitwise histogram control','required':True,
                       'passed':all(left['histograms_fb'][k]==right['histograms_fb'][k] for k in keys)})
    def compare(left,right,keys,label,required=True):
        values={}
        for key in keys:
            x,ex=left['equal_iteration_rates_pb'][key];y,ey=right['equal_iteration_rates_pb'][key]
            err=math.hypot(ex,ey)
            values[key]={'shift_pb':x-y,'combined_error_pb':err,'pull':(x-y)/err if err else (0 if x==y else float('inf'))}
        row={'label':label,'values':values,'passed':all(abs(v['pull'])<=3 for v in values.values()),'required':required,
             'kind':'common iteration estimates, separate-run combined errors; no cross-run event covariance'}
        checks.append(row)
        if not row['passed']:flags.append(row)
    for proc in a.processes:
        tags=['nominal','atomic','zero','grid8192','alpha01','slivers','seed2089','atomic2089']
        if (a.results/f'nlo_{proc}_cells_grid16384').exists():tags.append('grid16384')
        paths={tag:a.results/f'nlo_{proc}_cells_{tag}' for tag in tags}
        meta={tag:json.loads((path/'manifest.json').read_text()) for tag,path in paths.items()}
        for tag,m in meta.items():
            if m['status']!='completed' or m['exit_code']!=0:raise RuntimeError('Incomplete run: '+str(paths[tag]))
            if hashlib.sha256((paths[tag]/'input.ini').read_bytes()).hexdigest()!=m['config_sha256']:raise RuntimeError('Config hash mismatch')
            if m['executable_sha256']!=m['adapter_build']['adapter_executable_sha256']:raise RuntimeError('Executable mismatch')
            if m['environment'].get('LEPTONFF_FRACTION_CELLS')!=('0' if tag.startswith('atomic') else '1'):raise RuntimeError('Wrong fraction measure')
            if m['integration_warnings']:flags.append({'label':str(paths[tag]),'integration_warnings':m['integration_warnings']})
        raw={tag:analyse(path) for tag,path in paths.items()};raws[proc]=raw
        oldpath=reference/f'nlo_{proc}_coupled_nominal';oldmeta=json.loads((oldpath/'manifest.json').read_text())
        exact(meta['atomic'],oldmeta,list(oldmeta['histograms_fb']),proc+' inactive cells reproduce all v0.18 rates')
        exact(meta['nominal'],meta['atomic'],BASE,proc+' cells preserve continuum NLO and PHOTOS')
        exact(meta['grid8192'],meta['nominal'],BASE,proc+' grid leaves continuum NLO and PHOTOS invariant')
        if 'grid16384' in meta:exact(meta['grid16384'],meta['nominal'],BASE,proc+' finest grid leaves continuum NLO and PHOTOS invariant')
        exact(meta['seed2089'],meta['atomic2089'],BASE,proc+' independent seed continuum NLO and PHOTOS invariant')
        exact(meta['slivers'],meta['nominal'],['ff_born','ff_expanded','ff_matched','ff_order2','hard_projected','coupled_projected','photos_post'],
              proc+' geometric cutoffs leave FF/projected measurements invariant')
        zero=meta['zero']['histograms_fb'];born=zero['ff_born']
        passed=all(zero[k]==born for k in ['ff_expanded','ff_matched','decay_order1','coupled_matched'])
        passed=passed and all(zero[k][0][2:]==[0,0] for k in ['coupled_shift','coupled_order2','coupled_minus_order1'])
        passed=passed and zero['photos_post']==zero['photos_born'] and zero['photos_minus_coupled']==zero['photos_minus_ff']
        checks.append({'label':proc+' zero coupling, separate Born maps retained','kind':'exact','passed':passed,'required':True})
        compare(raw['alpha01'],raw['nominal'],['coupled_matched','coupled_shift','coupled_minus_order1','photos_minus_coupled'],proc+' QCD subtraction parameter')
        compare(raw['slivers'],raw['nominal'],['coupled_matched','coupled_shift','coupled_minus_order1','photos_minus_coupled'],proc+' geometric refinement')
        compare(raw['alpha01'],raw['nominal'],['coupled_projected','coupled_geometry','coupled_order2'],proc+' component diagnostic',False)
        compare(raw['nominal'],raw['atomic'],['coupled_minus_order1','coupled_order2'],proc+' finite-grid measure comparison',False)
        compare(raw['seed2089'],raw['nominal'],KEYS,proc+' independent-seed consistency',False)
        efficiency={}
        for cell,atom in [('nominal','atomic'),('seed2089','atomic2089')]:
            c=raw[cell]['equal_iteration_rates_pb'];r=raw[atom]['equal_iteration_rates_pb']
            tc,ta=cpu(paths[cell]),cpu(paths[atom]);stats={}
            for key in ['ff_beyond_order1','coupled_minus_order1','coupled_order2']:
                gain=(r[key][1]/c[key][1])**2
                stats[key]={'cell_pb':c[key],'atomic_pb':r[key],'observed_variance_gain_atomic_over_cell':gain,
                            'observed_cpu_efficiency_gain':gain*ta/tc}
            efficiency[cell]={'rates':stats,'cell_cpu_seconds':tc,'atomic_cpu_seconds':ta,
                              'paired_change':paired_chunks(paths[cell],paths[atom],KEYS)}
        refinement={tag:paired_chunks(paths[tag],paths['nominal'],KEYS) for tag in ['grid8192','slivers']+(['grid16384'] if 'grid16384' in paths else [])}
        oldgrid=paired_chunks(reference/f'nlo_{proc}_coupled_grid8192',oldpath,KEYS)
        summary[proc]={'nominal_rates_pb':{k:raw['nominal']['equal_iteration_rates_pb'][k] for k in KEYS},
                       'seed_efficiency_diagnostics':efficiency,'paired_cell_refinements':refinement,
                       'previous_atomic_grid8192_diagnostic':oldgrid}
        provenance[proc]={tag:{k:m[k] for k in ['source_sha256','working_tree_source_sha256','executable_sha256','config_sha256','environment','integration_warnings']}
                          for tag,m in meta.items()}
    failures=[c for c in checks if c['required'] and not c['passed']]
    report={'scope':'Conservative resolved-cell measurement on unchanged one-family NLL W/Z hard/FF kernels; endpoint atom exact',
            'checks':checks,'required_failures':failures,'diagnostic_flags':flags,'summary':summary,
            'raw_iteration_audits':raws,'provenance':provenance,
            'uncertainty_scope':'Observed variance/CPU ratios and paired ten-chunk changes are limited sample diagnostics; finite-grid target changes beyond order alpha; no universal variance or numerical error bound'}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'checks_passed':sum(c['passed'] for c in checks),'required_failures':len(failures),'diagnostic_flags':len(flags),'summary':summary},indent=2))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
