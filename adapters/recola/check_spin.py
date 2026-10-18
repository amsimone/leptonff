#!/usr/bin/env python3
"""Build and run the optional RECOLA 1.5.1 polarised tree-current check."""
import argparse,csv,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('recola',type=Path,help='Source installation with include/ and librecola.a')
p.add_argument('collier',type=Path,help='Directory containing libcollier.a')
p.add_argument('--core-build',type=Path,default=Path('build-host'))
p.add_argument('--cxx',default='/opt/homebrew/bin/g++-16')
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=Path(__file__).resolve().parents[2];out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
source=root/'adapters/recola/check_spin.cpp';core=a.core_build.resolve()/'libleptonff.a'
libraries=[a.recola.resolve()/'librecola.a',a.collier.resolve()/'libcollier.a']
exe=out/'check_spin'
cmd=[a.cxx,'-std=c++17','-O2','-DGNU_NAME_MANGLING',f'-I{root}/include',f'-I{a.recola.resolve()}/include',str(source),str(core),*[str(x) for x in libraries],'-lgfortran','-o',str(exe)]
with (out/'build.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
with (out/'run.log').open('w') as log:subprocess.run([str(exe),'points.csv'],cwd=out,stdout=log,stderr=subprocess.STDOUT,check=True)
rows=list(csv.DictReader((out/'points.csv').open()));maximum={}
for kind in ['0','1','2']:
    selected=[r for r in rows if r['kind']==kind]
    if len(selected)!=162:raise RuntimeError('Incomplete polarised grid')
    maximum[{'0':'wminus','1':'wplus','2':'neutral_parity_even'}[kind]]=max(abs(float(r['relative_difference'])) for r in selected)
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
report={'scope':'Massless real tree currents; both W charges and charge-symmetrised parity-even Z, three helicities, rotated and boosted kinematics',
        'points':len(rows),'maximum_relative_difference':maximum,'alpha':1/137.035999084,
        'widths':'zero, external stable bosons','normalisation':'unpolarised massless Born; M^2*(1-z)/(4*e^2) phase-space conversion',
        'command':cmd,'sha256':{str(path):sha(path) for path in [source,core,*libraries,exe,out/'points.csv']}}
(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if max(maximum.values())>1e-7:raise RuntimeError('Independent polarised current check failed')
