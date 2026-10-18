#!/usr/bin/env python3
"""Complete reproducible differential W/Z control and scale campaigns."""
import argparse,concurrent.futures,hashlib,json,subprocess,sys
from pathlib import Path

OPTIONS={'nominal':[],'zero':['--alpha','0'],'grid8192':['--bins','8192'],
         'scale_half':['--mu-frag','HALF'],'scale_double':['--mu-frag','DOUBLE'],
         'slivers':['--hard-soft','1e-6','--hard-collinear','1e-6'],'seed2089':['--seed','2089'],
         'qcd_half':['--mu-r-factor','.5','--mu-f-factor','.5'],
         'qcd_double':['--mu-r-factor','2','--mu-f-factor','2'],
         'qcd_rhalf':['--mu-r-factor','.5'],'qcd_rdouble':['--mu-r-factor','2'],
         'qcd_fhalf':['--mu-f-factor','.5'],'qcd_fdouble':['--mu-f-factor','2']}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('host',type=Path)
    p.add_argument('--adapter',type=Path,required=True);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--processes',nargs='+',choices=['wplus','wminus','z'],default=['wplus','wminus','z'])
    p.add_argument('--cases',nargs='+',choices=list(OPTIONS),default=list(OPTIONS))
    p.add_argument('--model',choices=['coupled','photonic'],default='coupled')
    p.add_argument('--call-scale',type=float,default=.1,help='Validation call budget; increase for precision')
    p.add_argument('--jobs',type=int,default=3);a=p.parse_args();a.results.mkdir(parents=True,exist_ok=True)
    if a.jobs<1:p.error('--jobs must be positive')
    runner=Path(__file__).with_name('run_reference.py')
    def run(process,case):
        out=a.results/f'nlo_{process}_spectra_{case}'
        if out.exists():
            meta=json.loads((out/'manifest.json').read_text())
            environment=meta['environment']
            if environment.get('LEPTONFF_PHOTONIC_LL')!=('1' if a.model=='photonic' else '0'):
                raise RuntimeError('Existing run uses a different theory model: '+str(out))
            if hashlib.sha256((out/'input.ini').read_bytes()).hexdigest()!=meta['config_sha256']:
                raise RuntimeError('Existing run configuration has changed: '+str(out))
            if meta['executable_sha256']!=hashlib.sha256(a.adapter.read_bytes()).hexdigest():
                raise RuntimeError('Existing run uses a different executable: '+str(out))
            if meta['status']!='completed':raise RuntimeError('Existing incomplete run: '+str(out))
            if meta.get('requested_inputs',{}).get('call_scale')!=a.call_scale:
                raise RuntimeError('Existing run has a different call budget: '+str(out))
            return process,case,'already completed'
        scale=91.1876 if process=='z' else 80.385
        options=[str(scale/2) if x=='HALF' else str(2*scale) if x=='DOUBLE' else x for x in OPTIONS[case]]
        command=[sys.executable,str(runner),str(a.host.resolve()),process,'--adapter',str(a.adapter.resolve()),
                 '--output',str(out.resolve()),'--threads','1','--spectra',
                 '--hard-samples','4','--call-scale',str(a.call_scale),
                 *(['--coupled-hard','--fraction-cells'] if a.model=='coupled' else ['--photonic-ll','--dress-hard']),*options]
        with Path(str(out)+'.launch.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        return process,case,'completed'
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        tasks=[pool.submit(run,process,case) for case in a.cases for process in a.processes]
        for task in concurrent.futures.as_completed(tasks):print(*task.result(),flush=True)
if __name__=='__main__':main()
