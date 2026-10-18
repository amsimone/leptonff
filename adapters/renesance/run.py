#!/usr/bin/env python3
"""Run isolated Born-QCD FSR references, preserving inputs and provenance."""
import argparse
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--client', type=Path, required=True)
    p.add_argument('--pdf-data', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--processes', nargs='+', choices=['z','wplus','wminus'], default=['z','wplus','wminus'])
    p.add_argument('--modes', nargs='+', choices=['fo','showers','analytic','all'], default=['fo','showers','analytic'])
    p.add_argument('--trials', type=int, default=1000000)
    p.add_argument('--analytic-trials', type=int, default=100000)
    p.add_argument('--seeds', type=int, nargs='+', default=[601,602,603,604])
    p.add_argument('--omega', type=float, default=1e-5)
    p.add_argument('--hard-samples', type=int, default=4)
    p.add_argument('--ff-bins', type=int, default=262144)
    p.add_argument('--frag-scale-factors', type=float, nargs='+', default=[1.0],
                   help='QED fragmentation scale / pole mass; scans require analytic mode alone')
    p.add_argument('--mass-matching',choices=['single_insertion','quadratic_finite'],default='single_insertion')
    p.add_argument('--parent-max', type=float, default=200)
    p.add_argument('--pythia-mec', choices=['true','false'], default='true')
    p.add_argument('--shower-cutoff', type=float, default=1e-5)
    p.add_argument('--mass-nodes', type=int, default=16)
    p.add_argument('--rapidity-nodes', type=int, default=64)
    p.add_argument('--real-nodes', type=int, default=8)
    p.add_argument('--z-monte-carlo', action='store_true', help='Use direct MC as an independent quadrature control')
    p.add_argument('--jobs', type=int, default=3)
    a = p.parse_args()
    if a.trials<2 or a.analytic_trials<2 or len(set(a.seeds))!=len(a.seeds) or any(s<0 or s>=1000000 for s in a.seeds):
        p.error('Need positive statistics and distinct independent seeds')
    if 'all' in a.modes and len(a.modes)!=1:
        p.error('all must be the only requested mode')
    if (len(set(a.frag_scale_factors))!=len(a.frag_scale_factors) or
            any(not math.isfinite(x) or x<=0 for x in a.frag_scale_factors)):
        p.error('Fragmentation scale factors must be finite, positive and distinct')
    if a.frag_scale_factors!=[1.0] and a.modes!=['analytic']:
        p.error('Scale scans require --modes analytic (avoid duplicated reference/shower runs)')
    if a.jobs<1:
        p.error('jobs must be positive')
    if a.mass_matching!='single_insertion' and (a.processes!=['z'] or a.modes!=['analytic']):
        p.error('Mass matching variations require --processes z --modes analytic')
    a.output.mkdir(parents=True, exist_ok=False)
    # A build during a campaign must not change later seeds' executable.
    client=a.output/'client'
    shutil.copy2(a.client,client)
    backend=a.client.resolve().parent/'manifest.json'
    if backend.exists():shutil.copy2(backend,a.output/'backend-manifest.json')
    env = dict(os.environ,LHAPDF_DATA_PATH=str(a.pdf_data.resolve()),OMP_NUM_THREADS='1')
    jobs=[]
    for process in a.processes:
        for mode in a.modes:
            for seed,frag_factor in ((s,f) for s in a.seeds for f in a.frag_scale_factors):
                zquad=process=='z' and mode in ('analytic','fo') and not a.z_monte_carlo
                if zquad and seed!=a.seeds[0]:
                    continue  # Deterministic quadrature is not independent seed statistics.
                name=f'{process}-{mode}-{seed}'
                if a.frag_scale_factors!=[1.0]:
                    name+= '-muf'+format(frag_factor,'.17g').replace('.','p')
                # Independent streams across separately integrated methods.
                # The on/off channels within a method still share each point.
                actual_seed=seed+1000000*(['fo','showers','analytic','all'].index(mode)+4*['z','wplus','wminus'].index(process))
                result=a.output/(name+'.json')
                cmd=[str(client.resolve()),str(result.resolve()),str({'z':202,'wplus':222,'wminus':212}[process]),
                     str(a.analytic_trials if mode=='analytic' else a.trials),str(actual_seed),str(a.omega),mode,
                     'NNPDF31_nnlo_as_0118',f'hard_samples={a.hard_samples}',f'ff_bins={a.ff_bins}',
                     f'm34max={a.parent_max}',f'pythia_mec={a.pythia_mec}',f'shower_cutoff={a.shower_cutoff}']
                if mode in ('analytic','all'):
                    cmd += [f'frag_scale_factor={frag_factor:.17g}']
                    if process=='z':cmd += [f'mass_matching={a.mass_matching}']
                if zquad:
                    cmd += ['z_quadrature=true',f'mass_nodes={a.mass_nodes}',f'rapidity_nodes={a.rapidity_nodes}',f'real_nodes={a.real_nodes}']
                jobs.append((name,cmd,result))
    manifest={'schema':'leptonff-born-fsr-campaign-v1','client_sha256':sha256(client),
              'pdf_data':str(a.pdf_data.resolve()),'argv':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
              'runs':[{'name':n,'command':c} for n,c,_ in jobs]}
    (a.output/'campaign.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def run(job):
        name,cmd,result=job
        with (a.output/(name+'.log')).open('w') as log:
            r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        if r.returncode:
            return name,False
        data=json.loads(result.read_text())
        data['provenance']={'client_sha256':manifest['client_sha256'],'command':cmd,
                            'pdf_data':manifest['pdf_data']}
        result.write_text(json.dumps(data,separators=(',',':'))+'\n')
        return name,True
    failed=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for name,ok in pool.map(run,jobs):
            print(name,'complete' if ok else 'FAILED',flush=True)
            if not ok:failed.append(name)
    if failed:
        raise SystemExit('Failed runs: '+', '.join(failed))


if __name__=='__main__':
    main()
