#!/usr/bin/env python3
"""Run immutable all-parent coupled hard controls with raw iteration output."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('host',type=Path);p.add_argument('process',choices=['wplus','wminus','z'])
p.add_argument('--adapter',type=Path,required=True);p.add_argument('--results',type=Path,default=Path('results'))
p.add_argument('--cases',nargs='+',choices=['nominal','zero','grid8192','grid16384','alpha01','slivers','scale_half','scale_double'],
               default=['nominal','zero','grid8192','alpha01','slivers','scale_half','scale_double'])
a=p.parse_args();root=Path(__file__).resolve().parents[2]
M=91.1876 if a.process=='z' else 80.385
options={'nominal':[],'zero':['--alpha','0'],'grid8192':['--bins','8192'],'grid16384':['--bins','16384'],'alpha01':['--dipole-alpha','.1'],
         'slivers':['--hard-soft','1e-6','--hard-collinear','1e-6'],
         'scale_half':['--mu-frag',str(M/2)],'scale_double':['--mu-frag',str(M*2)]}
for case in a.cases:
    out=a.results/f'nlo_{a.process}_coupled_{case}'
    command=[sys.executable,str(root/'adapters/mcfm/run_reference.py'),str(a.host.resolve()),a.process,
             '--output',str(out.resolve()),'--adapter',str(a.adapter.resolve()),'--threads','1',
             '--hard-samples','4','--coupled-hard',*options[case]]
    with Path(str(out)+'.launch.log').open('w') as log:
        subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    print(case,'completed',flush=True)
