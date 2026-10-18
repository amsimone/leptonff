#!/usr/bin/env python3
"""Paired first-order PHOTOS proposal test, physical alpha versus alpha/4."""
import argparse,hashlib,json,math,os,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('executable',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--events',type=int,default=50000)
a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=False);results=[]
for proc in [1,6,31]:
    branches=[]
    for scale in [1,4]:
        path=out/f'process{proc}_scale{scale}.txt'
        env={**os.environ,'OMP_NUM_THREADS':'1','LEPTONFF_PHOTOS_ORDER':'1',
             'LEPTONFF_PHOTOS_CUTOFF':'1e-3','LEPTONFF_PHOTOS_LINEAR_SCALE':str(scale),
             'LEPTONFF_ALPHA':str(1/137.035999084)}
        with (out/f'process{proc}_scale{scale}.log').open('w') as log:
            subprocess.run([str(a.executable.resolve()),str(proc),str(a.events),str(path.resolve())],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        rows=[float(v) for v in path.read_text().splitlines()]
        if len(rows)!=a.events:raise RuntimeError('Incomplete proposal output')
        branches.append(rows)
    differences=[x-y for x,y in zip(*branches)];n=len(differences);mean=sum(differences)/n
    se=math.sqrt(sum((x-mean)**2 for x in differences)/(n*(n-1)))
    pull=mean/se if se else 0
    results.append({'process':proc,'events':n,'scale1_mean':sum(branches[0])/n,
                    'scale4_mean':sum(branches[1])/n,'paired_difference':mean,'paired_error':se,'pull':pull})
report={'checks':results,'executable_sha256':hashlib.sha256(a.executable.read_bytes()).hexdigest(),
        'scope':'Statistical linear-proposal closure at rest; not a matrix-element accuracy comparison'}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if any(abs(r['pull'])>5 for r in results):raise RuntimeError('PHOTOS linear proposal check failed')
