#!/usr/bin/env python3
"""Audit the one-family coupled W/Z hard pilot and retain rate diagnostics."""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path
from analyse_iterations import analyse


def paired_chunks(left,right,keys):
    """Paired iteration-mean diagnostic; event cross-products are unavailable."""
    def rows(path):
        out=defaultdict(dict)
        with (path/'iteration_rates.csv').open() as file:
            for row in csv.reader(file):
                stage,part,ips,it,calls=map(int,row[:5])
                if stage!=1:raise RuntimeError('Unexpected stage')
                out[(part,ips,it)][row[5]]=(float(row[6]),calls)
        return out
    a,b=rows(left),rows(right)
    if a.keys()!=b.keys():raise RuntimeError('Cannot pair different iteration sets')
    grouped=defaultdict(list)
    for index in a:
        if a[index]['ff_born']!=b[index]['ff_born']:
            raise RuntimeError('Paired refinement requires identical Born iteration control')
        if any(a[index][k][1]!=b[index][k][1] for k in keys):raise RuntimeError('Paired calls differ')
        grouped[index[:2]].append({k:(a[index][k][0]-b[index][k][0])*1e-3 for k in keys})
    out={k:{'shift_pb':0.,'paired_chunk_error_pb':0.} for k in keys}
    for data in grouped.values():
        n=len(data)
        if n<2:raise RuntimeError('Too few paired chunks')
        for key in keys:
            mean=sum(d[key] for d in data)/n
            out[key]['shift_pb']+=mean
            out[key]['paired_chunk_error_pb']+=sum((d[key]-mean)**2 for d in data)/(n*(n-1))
    for value in out.values():value['paired_chunk_error_pb']=math.sqrt(value['paired_chunk_error_pb'])
    return {'values':out,'chunks_per_contribution':{str(k):len(v) for k,v in grouped.items()},
            'scope':'Paired equal-call raw iteration means, identical Born controls; limited chunk-error diagnostic, no event-level cross-products or heavy-tail bound'}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,default=Path('results'))
    p.add_argument('--output',type=Path,default=Path('results/coupled_lab_host.json'))
    p.add_argument('--processes',nargs='+',choices=['wplus','wminus','z'],default=['wplus','wminus','z'])
    a=p.parse_args();checks=[];flags=[];rates={};raws={};provenance={}
    def exact(left,right,keys,label):
        passed=all(left['histograms_fb'][key]==right['histograms_fb'][key] for key in keys)
        checks.append({'label':label,'kind':'same-stream exact','passed':passed,'required':True})
    def comparison(left,right,keys,label,required=True):
        values={}
        for key in keys:
            x,ex=left['equal_iteration_rates_pb'][key];y,ey=right['equal_iteration_rates_pb'][key]
            error=math.hypot(ex,ey);pull=(x-y)/error if error else (0 if x==y else float('inf'))
            values[key]={'shift_pb':x-y,'combined_error_pb':error,'pull':pull}
        row={'label':label,'kind':'common iteration means; separate-run combined errors, no cross-run covariance',
             'values':values,'passed':all(abs(v['pull'])<=3 for v in values.values()),'required':required}
        checks.append(row)
        if not row['passed']:flags.append(row)
    for proc in a.processes:
        paths={tag:a.results/f'nlo_{proc}_coupled_{tag}' for tag in
               ('nominal','zero','grid8192','alpha01','slivers','scale_half','scale_double')}
        if (a.results/f'nlo_{proc}_coupled_grid16384/manifest.json').exists():
            paths['grid16384']=a.results/f'nlo_{proc}_coupled_grid16384'
        meta={tag:json.loads((path/'manifest.json').read_text()) for tag,path in paths.items()}
        for tag,m in meta.items():
            if m['executable_sha256']!=m['adapter_build']['adapter_executable_sha256']:raise RuntimeError('Executable hash mismatch')
            if m['environment'].get('LEPTONFF_COUPLED_HARD')!='1':raise RuntimeError('Wrong hard matching mode')
            if m['status']!='completed' or m['exit_code']!=0:raise RuntimeError('Incomplete run: '+str(paths[tag]))
            if hashlib.sha256((paths[tag]/'input.ini').read_bytes()).hexdigest()!=m['config_sha256']:raise RuntimeError('Config hash differs')
            if m['integration_warnings']:flags.append({'label':str(paths[tag]),'integration_warnings':m['integration_warnings']})
        raw={tag:analyse(path) for tag,path in paths.items()};raws[proc]=raw
        nominal=meta['nominal'];zero=meta['zero'];grid=meta['grid8192'];fine=meta['slivers']
        old=json.loads((a.results/f'nlo_{proc}_dressed_audit/manifest.json').read_text())
        baseline=['ff_born','ff_expanded','hard_projected','hard_geometry','hard_shift','decay_order1',
                  'photos_born','photos_post','photos_shift','photos_radiated','photos_minus_order1']
        exact(nominal,old,baseline,proc+' preserves continuum NLO and PHOTOS baseline')
        ff_baseline=a.results/('nlo_z_coupled_ff_control' if proc=='z' else f'nlo_{proc}_photos_seeded')
        oldff=json.loads((ff_baseline/'manifest.json').read_text())
        exact(nominal,oldff,['ff_born','ff_expanded','ff_resummed','ff_matched','ff_shift','ff_beyond_order1','ff_order2'],
              proc+' preserves coupled FF-only baseline')
        identities=['ff_expanded','ff_matched','decay_order1','coupled_matched']
        passed=all(zero['histograms_fb'][key]==zero['histograms_fb']['ff_born'] for key in identities)
        passed=passed and all(zero['histograms_fb'][key][0][2:]==[0,0] for key in
                             ['coupled_projected','coupled_geometry','coupled_shift','coupled_order2','coupled_minus_order1'])
        passed=passed and zero['histograms_fb']['photos_post']==zero['histograms_fb']['photos_born']
        passed=passed and zero['histograms_fb']['photos_minus_coupled']==zero['histograms_fb']['photos_minus_ff']
        checks.append({'label':proc+' zero-coupling identities','kind':'exact','passed':passed,'required':True})
        exact(nominal,grid,baseline,proc+' FF grid leaves continuum NLO and PHOTOS unchanged')
        if 'grid16384' in meta:
            exact(nominal,meta['grid16384'],baseline,proc+' 16384 grid leaves continuum NLO and PHOTOS unchanged')
        exact(nominal,fine,['ff_born','ff_expanded','ff_matched','hard_projected','coupled_projected','photos_post'],
              proc+' geometric slivers leave FF/projected sources unchanged')
        comparison(raw['alpha01'],raw['nominal'],['coupled_matched','coupled_shift','coupled_minus_order1','photos_minus_coupled'],proc+' QCD subtraction parameter')
        comparison(raw['slivers'],raw['nominal'],['coupled_matched','coupled_shift','coupled_minus_order1','photos_minus_coupled'],proc+' finite geometric refinement')
        comparison(raw['alpha01'],raw['nominal'],['coupled_projected','coupled_geometry','coupled_order2'],proc+' component diagnostic',False)
        # Factorisation-scale changes cancel analytically at continuum NLO.
        # Printed MCFM rows can differ in their last rounding digit.
        scale_controls={}
        for tag in ('scale_half','scale_double'):
            values={}
            for key in ('decay_order1','ff_born','photos_post'):
                l=meta[tag]['histograms_fb'][key][0][2:];r=nominal['histograms_fb'][key][0][2:]
                residual=max(abs(x-y)*1e-3 for x,y in zip(l,r));values[key]=residual
            scale_controls[tag]=values
        checks.append({'label':proc+' continuum NLO fragmentation-scale cancellation','kind':'printed means/errors within 1e-6 pb',
                       'values':scale_controls,'passed':all(v<=1e-6 for row in scale_controls.values() for v in row.values()),'required':True})
        rates[proc]={key:raw['nominal']['equal_iteration_rates_pb'][key] for key in
                     ['coupled_projected','coupled_geometry','coupled_shift','coupled_matched','coupled_order2',
                      'coupled_minus_order1','photos_minus_coupled','decay_order1','decay_matched']}
        rates[proc]['refinement_shifts_pb']={tag:{key:raw[tag]['equal_iteration_rates_pb'][key][0]-raw['nominal']['equal_iteration_rates_pb'][key][0]
                  for key in ['coupled_matched','coupled_shift','coupled_order2','coupled_minus_order1']} for tag in
                  (['grid8192','slivers','scale_half','scale_double']+(['grid16384'] if 'grid16384' in raw else []))}
        rates[proc]['paired_refinement_diagnostics']={tag:paired_chunks(paths[tag],paths['nominal'],
            ['coupled_matched','coupled_minus_order1','coupled_order2']) for tag in
            (['grid8192','slivers','scale_half','scale_double']+(['grid16384'] if 'grid16384' in raw else []))}
        provenance[proc]={tag:{**{key:m[key] for key in ['source_sha256','executable_sha256','config_sha256','environment','integration_warnings']},
                          'build_source_sha256':m['adapter_build']['source_sha256'],
                          'launch_working_tree_source_sha256':m.get('working_tree_source_sha256',m['source_sha256'])}
                          for tag,m in meta.items()}
    failures=[c for c in checks if c['required'] and not c['passed']]
    report={'scope':'One-family MSbar NLL hard/FF spin W/Z on NLO QCD; physical flavours/thresholds absent; PHOTOS diagnostic only',
            'rates_pb':rates,'checks':checks,'diagnostic_flags':flags,'required_failures':failures,
            'raw_iteration_audits':raws,'provenance':provenance,
            'uncertainty_scope':'Common iteration means and within-iteration errors; ten-chunk covariance diagnostic; no certified physics, heavy-tail or cross-bin bound'}
    a.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'checks_passed':sum(c['passed'] for c in checks),'required_failures':len(failures),'diagnostic_flags':len(flags),'rates_pb':rates},indent=2))
    if failures:raise SystemExit(1)

if __name__=='__main__':main()
