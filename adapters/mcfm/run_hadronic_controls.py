#!/usr/bin/env python3
"""Sourced HVP/LO-QCD pointlike-boundary W/Z controls; separate legacy regression."""
import argparse,concurrent.futures,hashlib,json,subprocess,sys
from pathlib import Path

OPTIONS={'nominal':[],'zero':['--alpha','0'],'qcd_off':['--qcd-off'],
    'had_minus':['--hadronic-shift','-1'],'had_plus':['--hadronic-shift','1'],
    'match_low':['--hadronic-scale','1.5'],'match_high':['--hadronic-scale','3'],
    'grid1024':['--bins','1024'],
    'alternative':['--hadronic-input','ALTERNATIVE'],
    'scale_half':['--mu-frag','HALF'],'scale_double':['--mu-frag','DOUBLE'],
    'legacy':['--flavours','one-family','--bins','4096'],
    'leptons':['--flavours','leptons']}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('host',type=Path)
    p.add_argument('--adapter',type=Path,required=True);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--alternative',type=Path,default=Path('results/external_hvp_v022/jegerlehner19.dat'))
    p.add_argument('--cases',nargs='+',choices=list(OPTIONS),default=list(OPTIONS))
    p.add_argument('--processes',nargs='+',choices=['wplus','wminus','z'],default=['wplus','wminus','z'])
    p.add_argument('--jobs',type=int,default=3);p.add_argument('--call-scale',type=float,default=.05)
    a=p.parse_args();a.results.mkdir(parents=True,exist_ok=True)
    if a.jobs<1:p.error('jobs must be positive')
    exe_hash=hashlib.sha256(a.adapter.read_bytes()).hexdigest()
    def run(process,case):
        out=a.results/f'nlo_{process}_hadronic_{case}';M=91.1876 if process=='z' else 80.385
        opts=[str(M/2) if v=='HALF' else str(2*M) if v=='DOUBLE' else str(a.alternative.resolve()) if v=='ALTERNATIVE' else v for v in OPTIONS[case]]
        command=[sys.executable,str(Path(__file__).with_name('run_reference.py')),str(a.host.resolve()),process,
            '--adapter',str(a.adapter.resolve()),'--output',str(out.resolve()),'--threads','1',
            '--spectra','--coupled-hard','--fraction-cells','--hard-samples','4','--call-scale',str(a.call_scale),
            '--flavours','hadronic-pointlike','--bins','512',*opts]
        if out.exists():
            m=json.loads((out/'manifest.json').read_text())
            if m['status']!='completed' or m['executable_sha256']!=exe_hash or json.loads((out/'launch.json').read_text())!=command:
                raise RuntimeError('Existing incomplete or different run: '+str(out))
            if hashlib.sha256((out/'input.ini').read_bytes()).hexdigest()!=m['config_sha256']:raise RuntimeError('Changed configuration')
            if case not in ['legacy','leptons']:
                source=a.alternative if case=='alternative' else Path(__file__).resolve().parents[2]/'configs/hadronic_alpha_pdg2025.dat'
                if hashlib.sha256(source.read_bytes()).hexdigest()!=m['hadronic_input']['file_sha256'] or hashlib.sha256((out/'hadronic.dat').read_bytes()).hexdigest()!=m['hadronic_input']['file_sha256']:
                    raise RuntimeError('Changed hadronic input: '+str(out))
            return process,case,'already completed'
        with Path(str(out)+'.launch.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        (out/'launch.json').write_text(json.dumps(command,indent=2)+'\n')
        return process,case,'completed'
    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        futures=[pool.submit(run,process,case) for case in a.cases for process in a.processes]
        for f in concurrent.futures.as_completed(futures):print(*f.result(),flush=True)
if __name__=='__main__':main()
