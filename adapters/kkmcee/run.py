#!/usr/bin/env python3
"""Run isolated, pinned fixed-Q KKMC references with paired perturbative orders."""
import argparse
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess

HERE=Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--client',type=Path,default=HERE.parents[1]/'build-kkmcee-fsr/kkmcee_fsr_reference')
    p.add_argument('--defaults',type=Path,default=HERE.parents[1]/'external/kkmcee-f77-audit/.KK2f_defaults')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--Q',type=float,nargs='+',default=[75,91.1876,120])
    p.add_argument('--mass',type=float,default=.1056583755)
    p.add_argument('--alpha-factor',type=float,default=1)
    p.add_argument('--trials',type=int,default=1000000)
    p.add_argument('--seed',type=int,default=913001)
    p.add_argument('--seeds',type=int,nargs='+',help='Independent base seeds; each Q adds its index times 100003')
    p.add_argument('--cutoff',type=float,default=1e-5)
    p.add_argument('--fsr-off',action='store_true')
    p.add_argument('--jobs',type=int,default=1)
    a=p.parse_args()
    seeds=a.seeds if a.seeds is not None else [a.seed]
    resolved=[seed+index*100003 for seed in seeds for index in range(len(a.Q))]
    if (a.trials<2 or a.jobs<1 or len(set(a.Q))!=len(a.Q) or
            len(set(resolved))!=len(resolved) or any(s<=0 or s>=900000000 for s in resolved) or
            any(not math.isfinite(q) or q<=10*a.mass for q in a.Q) or
            not math.isfinite(a.alpha_factor) or not 0<a.alpha_factor<10):
        p.error('Invalid statistics, seeds, coupling or parent masses')
    a.output.mkdir(parents=True,exist_ok=False)
    executable=a.output/'client';shutil.copy2(a.client,executable)
    defaults=a.output/'defaults.input';shutil.copy2(a.defaults,defaults)
    backend=a.client.resolve().parent/'manifest.json'
    if not backend.exists():backend=a.client.resolve().parent/'backend-manifest.json'
    if backend.exists():shutil.copy2(backend,a.output/'backend-manifest.json')
    manifest=dict(schema='leptonff-kkmcee-campaign-v1',client_sha256=digest(executable),
                  defaults_sha256=digest(defaults),runs=[])
    for seed in seeds:
        for index,Q in enumerate(a.Q):
            name='Q'+format(Q,'.8g').replace('.','p')
            if len(seeds)>1:name+='-seed'+str(seed+index*100003)
            folder=a.output/name;folder.mkdir()
            result=folder/'result.json'
            cmd=[str(executable.resolve()),str(defaults.resolve()),str(result.resolve()),
                 str(Q),str(a.mass),str(a.alpha_factor/137.035999084),str(a.trials),
                 str(seed+index*100003),str(a.cutoff),'0' if a.fsr_off else '1']
            manifest['runs'].append(dict(name=name,command=cmd))
    (a.output/'campaign.json').write_text(json.dumps(manifest,indent=2)+'\n')
    def run(job):
        folder=a.output/job['name'];result=folder/'result.json'
        with (folder/'run.log').open('w') as log:
            r=subprocess.run(job['command'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,
                             env=dict(os.environ,OMP_NUM_THREADS='1'))
        # Legacy Fortran STOP can return zero: the complete result is mandatory.
        if r.returncode or not result.exists():
            raise RuntimeError('KKMC failed: '+job['name']+'; inspect run.log/fort.16')
        d=json.loads(result.read_text())
        if d['schema']!='leptonff-kkmcee-fixed-q-v1' or d['trials']!=a.trials:
            raise ValueError('Invalid KKMC result')
        d['provenance']=dict(client_sha256=manifest['client_sha256'],
                             defaults_sha256=manifest['defaults_sha256'],command=job['command'])
        result.write_text(json.dumps(d,separators=(',',':'))+'\n')
        return job['name'],d['seconds']
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for name,seconds in pool.map(run,manifest['runs']):
            print(name,'complete',round(seconds,2),'seconds',flush=True)


if __name__=='__main__':main()
