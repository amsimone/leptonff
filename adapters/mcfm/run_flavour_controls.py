#!/usr/bin/env python3
"""Reproducible variable-flavour W/Z matching controls, including legacy mode."""
import argparse,concurrent.futures,hashlib,json,subprocess,sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
OPTIONS={
    'nominal':[], 'threshold_half':['--threshold-factor','.5'],
    'threshold_double':['--threshold-factor','2'], 'zero':['--alpha','0'],
    'msbar':['--alpha-scheme','msbar'], 'legacy':['--flavours','one-family'],
    'file':['--flavours','custom','--flavour-file',str(ROOT/'configs/leptons_thresholds.dat'),
            '--alpha-scheme','on-shell'],
    'grid8192':['--bins','8192'], 'scale_half':['--mu-frag','HALF'],
    'scale_double':['--mu-frag','DOUBLE'],
    'slivers':['--hard-soft','1e-6','--hard-collinear','1e-6'],
    'test_charges':['--flavours','custom','--flavour-file',str(ROOT/'configs/charged_fermion_threshold_test.dat'),
                    '--alpha-scheme','msbar']}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('host',type=Path)
    p.add_argument('--adapter',type=Path,required=True);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--processes',nargs='+',choices=['wplus','wminus','z'],default=['wplus','wminus','z'])
    p.add_argument('--cases',nargs='+',choices=list(OPTIONS),default=list(OPTIONS))
    p.add_argument('--call-scale',type=float,default=.05);p.add_argument('--jobs',type=int,default=3)
    a=p.parse_args();a.results.mkdir(parents=True,exist_ok=True)
    if a.jobs<1:p.error('--jobs must be positive')
    exe_hash=hashlib.sha256(a.adapter.read_bytes()).hexdigest()
    def run(process,case):
        out=a.results/f'nlo_{process}_flavours_{case}'
        scale=91.1876 if process=='z' else 80.385
        opts=[str(scale/2) if x=='HALF' else str(2*scale) if x=='DOUBLE' else x for x in OPTIONS[case]]
        command=[sys.executable,str(Path(__file__).with_name('run_reference.py')),str(a.host.resolve()),process,
                 '--adapter',str(a.adapter.resolve()),'--output',str(out.resolve()),'--threads','1',
                 '--spectra','--coupled-hard','--fraction-cells','--hard-samples','4',
                 '--call-scale',str(a.call_scale),'--flavours','leptons',*opts]
        if out.exists():
            m=json.loads((out/'manifest.json').read_text())
            if m['status']!='completed' or m['executable_sha256']!=exe_hash:
                raise RuntimeError('Incomplete or different executable: '+str(out))
            if m['requested_inputs']['call_scale']!=a.call_scale:
                raise RuntimeError('Different call budget: '+str(out))
            if hashlib.sha256((out/'input.ini').read_bytes()).hexdigest()!=m['config_sha256']:
                raise RuntimeError('Modified configuration: '+str(out))
            if json.loads((out/'launch.json').read_text())!=command:
                raise RuntimeError('Different launch options: '+str(out))
            return process,case,'already completed'
        with Path(str(out)+'.launch.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        (out/'launch.json').write_text(json.dumps(command,indent=2)+'\n')
        return process,case,'completed'
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        tasks=[pool.submit(run,process,case) for case in a.cases for process in a.processes]
        for task in concurrent.futures.as_completed(tasks):print(*task.result(),flush=True)

if __name__=='__main__':main()
