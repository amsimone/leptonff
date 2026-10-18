#!/usr/bin/env python3
"""Run signed conservative-cell W/Z host refinement controls."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('host',type=Path);p.add_argument('process',choices=['wplus','wminus','z'])
p.add_argument('--adapter',type=Path,required=True);p.add_argument('--results',type=Path,default=Path('results'))
p.add_argument('--cases',nargs='+',choices=['nominal','atomic','zero','grid8192','grid16384','alpha01','slivers','seed2089','atomic2089'],
               default=['nominal','atomic','zero','grid8192','alpha01','slivers'])
a=p.parse_args();root=Path(__file__).resolve().parents[2]
a.results.mkdir(parents=True,exist_ok=True)
options={'nominal':[],'atomic':[],'zero':['--alpha','0'],'grid8192':['--bins','8192'],
         'grid16384':['--bins','16384'],'alpha01':['--dipole-alpha','.1'],
         'slivers':['--hard-soft','1e-6','--hard-collinear','1e-6'],'seed2089':['--seed','2089'],'atomic2089':['--seed','2089']}
for case in a.cases:
    out=a.results/f'nlo_{a.process}_cells_{case}'
    command=[sys.executable,str(root/'adapters/mcfm/run_reference.py'),str(a.host.resolve()),a.process,
             '--output',str(out.resolve()),'--adapter',str(a.adapter.resolve()),'--threads','1',
             '--hard-samples','4','--coupled-hard',*options[case]]
    if not case.startswith('atomic'):command.append('--fraction-cells')
    with Path(str(out)+'.launch.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(case,'completed',flush=True)
