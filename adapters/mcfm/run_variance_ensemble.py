#!/usr/bin/env python3
"""Independent W-minus seed ensemble at smaller QCD dipole alpha."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from run_comparison import MODELS

SEEDS = [1045, 2089, 3791, 5003, 7207, 9001, 11003, 13007, 15013, 17011]


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('host', type=Path)
    p.add_argument('--adapter', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=6)
    p.add_argument('--call-scale', type=float, default=.5)
    a = p.parse_args()
    if a.jobs < 1:
        p.error('jobs must be positive')
    a.results.mkdir(parents=True, exist_ok=True)
    executable = digest(a.adapter)
    runs = {model: {'wminus': [str((a.results/f'{model}_wminus_seed{seed}').resolve())
                               for seed in SEEDS]} for model in MODELS}
    (a.results/'central_runs.json').write_text(json.dumps(runs, indent=2)+'\n')

    def run(model, seed):
        out = a.results/f'{model}_wminus_seed{seed}'
        bins, span = (1024, 64) if model == 'qed_hvp' else (4096, 40)
        command = [sys.executable, str(Path(__file__).with_name('run_reference.py').resolve()),
                   str(a.host.resolve()), 'wminus', '--adapter', str(a.adapter.resolve()),
                   '--output', str(out.resolve()), '--threads', '1', '--spectra',
                   '--hard-samples', '16', '--call-scale', str(a.call_scale),
                   '--bins', str(bins), '--log-span', str(span), '--seed', str(seed),
                   '--dipole-alpha', '.01', *MODELS[model]]
        if out.exists():
            meta = json.loads((out/'manifest.json').read_text())
            if (meta['status'] != 'completed' or meta['executable_sha256'] != executable
                    or json.loads((out/'launch.json').read_text()) != command
                    or digest(out/'input.ini') != meta['config_sha256']):
                raise RuntimeError('Existing incomplete or different run: '+str(out))
            if model == 'qed_hvp':
                source = Path(__file__).resolve().parents[2]/'configs/hadronic_alpha_pdg2025.dat'
                if digest(source) != meta['hadronic_input']['file_sha256'] or digest(out/'hadronic.dat') != digest(source):
                    raise RuntimeError('Changed hadronic input: '+str(out))
            return model, seed, 'already completed'
        with Path(str(out)+'.launch.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        (out/'launch.json').write_text(json.dumps(command, indent=2)+'\n')
        return model, seed, 'completed'

    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        futures = [pool.submit(run, model, seed) for seed in SEEDS for model in MODELS]
        for future in concurrent.futures.as_completed(futures):
            print(*future.result(), flush=True)


if __name__ == '__main__':
    main()
