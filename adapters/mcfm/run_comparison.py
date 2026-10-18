#!/usr/bin/env python3
"""Reproducible pure-QED W/Z comparison, with explicit numerical controls."""
import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import subprocess
import sys

MODELS = {
    'photonic': ['--photonic-ll', '--dress-hard'],
    'photonic_nll': ['--photonic-nll', '--dress-hard'],
    'leptons': ['--coupled-hard', '--fraction-cells', '--flavours', 'leptons'],
    # Existing v0.22 input prescription, with strong FF evolution disabled.
    # This remains an explicit pointlike hadronic boundary approximation.
    'qed_hvp': ['--coupled-hard', '--fraction-cells', '--flavours',
                'hadronic-pointlike', '--qcd-off'],
}
CASES = {
    'nominal': [], 'seed': ['--seed', '2089'],
    'grid': ['--bins', 'DOUBLE_GRID'], 'domain': ['--log-span', 'WIDE_DOMAIN'],
    'tail': ['--log-span', 'DOUBLE_DOMAIN', '--bins', 'DOUBLE_GRID'],
    'angular': ['--hard-samples', '16'],
    'angular_seed': ['--hard-samples', '16', '--seed', '2089'],
    'cutoff': ['--hard-soft', '1e-6', '--hard-collinear', '1e-6'],
    'frag_half': ['--mu-frag', 'HALF'], 'frag_double': ['--mu-frag', 'DOUBLE'],
    'photos_cutoff': ['--photos-cutoff', '1e-6'],
    'qcd_half': ['--mu-r-factor', '.5', '--mu-f-factor', '.5'],
    'qcd_double': ['--mu-r-factor', '2', '--mu-f-factor', '2'],
    'qcd_rhalf': ['--mu-r-factor', '.5'], 'qcd_rdouble': ['--mu-r-factor', '2'],
    'qcd_fhalf': ['--mu-f-factor', '.5'], 'qcd_fdouble': ['--mu-f-factor', '2'],
    'threshold_low': ['--hadronic-scale', '1.5'],
    'threshold_high': ['--hadronic-scale', '3'],
    'hvp_low': ['--hadronic-shift', '-1'], 'hvp_high': ['--hadronic-shift', '1'],
    'dipole_tenth': ['--dipole-alpha', '.1', '--hard-samples', '16'],
    'dipole_tenth_seed': ['--dipole-alpha', '.1', '--hard-samples', '16', '--seed', '2089'],
    'dipole_hundredth': ['--dipole-alpha', '.01', '--hard-samples', '16'],
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('host', type=Path)
    p.add_argument('--adapter', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--models', nargs='+', choices=MODELS, default=['photonic', 'leptons'])
    p.add_argument('--processes', nargs='+', choices=['wplus', 'wminus', 'z'],
                   default=['wplus', 'wminus', 'z'])
    p.add_argument('--cases', nargs='+', choices=CASES,
                   default=['nominal', 'seed', 'grid', 'domain', 'angular', 'cutoff'])
    p.add_argument('--call-scale', type=float, default=.5)
    p.add_argument('--jobs', type=int, default=3)
    p.add_argument('--pdf', default='NNPDF31_nnlo_as_0118')
    a = p.parse_args()
    if a.jobs < 1:
        p.error('--jobs must be positive')
    hadronic_cases = {'threshold_low', 'threshold_high', 'hvp_low', 'hvp_high'}
    if hadronic_cases.intersection(a.cases) and a.models != ['qed_hvp']:
        p.error('Hadronic input controls require --models qed_hvp alone')
    a.results.mkdir(parents=True, exist_ok=True)
    executable_hash = digest(a.adapter)

    def run(model, process, case):
        out = a.results / f'{model}_{process}_{case}'
        mass = 91.1876 if process == 'z' else 80.385
        bins, span = (1024, 64) if model == 'qed_hvp' else (4096, 40)
        substitutions = {'DOUBLE_GRID': str(2*bins), 'WIDE_DOMAIN': str(span*1.5),
                         'DOUBLE_DOMAIN': str(2*span),
                         'HALF': str(mass/2), 'DOUBLE': str(2*mass)}
        options = [substitutions.get(x, x) for x in CASES[case]]
        if model == 'qed_hvp' and case == 'domain':
            options = ['--log-span', '80']
        if model == 'qed_hvp' and case == 'tail':
            # A shorter, still bound-checked domain preserves h=span/bins.
            # Doubling 64 would exceed the coupled backend's 80 limit.
            options = ['--log-span', '32', '--bins', str(bins//2)]
        command = [sys.executable, str(Path(__file__).with_name('run_reference.py').resolve()),
                   str(a.host.resolve()), process, '--adapter', str(a.adapter.resolve()),
                   '--output', str(out.resolve()), '--threads', '1', '--spectra',
                   '--hard-samples', '4', '--call-scale', str(a.call_scale),
                   '--bins', str(bins), '--log-span', str(span), '--pdf', a.pdf,
                   *MODELS[model], *options]
        if out.exists():
            meta = json.loads((out/'manifest.json').read_text())
            if (meta['status'] != 'completed' or meta['executable_sha256'] != executable_hash
                    or json.loads((out/'launch.json').read_text()) != command
                    or digest(out/'input.ini') != meta['config_sha256']):
                raise RuntimeError('Incomplete or different existing run: '+str(out))
            if model == 'qed_hvp':
                source = Path(__file__).resolve().parents[2]/'configs/hadronic_alpha_pdg2025.dat'
                if digest(source) != meta['hadronic_input']['file_sha256'] or digest(out/'hadronic.dat') != digest(source):
                    raise RuntimeError('Changed hadronic input: '+str(out))
            return model, process, case, 'already completed'
        with Path(str(out)+'.launch.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        (out/'launch.json').write_text(json.dumps(command, indent=2)+'\n')
        return model, process, case, 'completed'

    with concurrent.futures.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        futures = [pool.submit(run, model, process, case)
                   for case in a.cases for model in a.models for process in a.processes]
        for future in concurrent.futures.as_completed(futures):
            print(*future.result(), flush=True)


if __name__ == '__main__':
    main()
