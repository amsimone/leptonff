#!/usr/bin/env python3
"""Run independently seeded McMule pieces, retaining signed shape integrals."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import time

PIECES = ('ee2mm0', 'ee2mmFMM', 'ee2mmRMM', 'ee2mmRFMMMM', 'ee2mmRRMMMM')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(build, manifest):
    for name, digest in manifest['artifacts'].items():
        if sha(build / name) != digest:
            raise ValueError('Changed build artifact: ' + name)
    for name, digest in manifest['runtime_artifacts'].items():
        if sha(name) != digest:
            raise ValueError('Changed Fortran runtime: ' + name)


def execute(config):
    folder = Path(config['output']); folder.mkdir(parents=True, exist_ok=False)
    build = Path(config['build']); manifest = json.loads((build / 'manifest.json').read_text())
    verify(build, manifest)
    args = [str(build / Path(manifest['executable']).name), 'result.json', config['piece'], config['mode']]
    args += [str(config[k]) for k in ('Q', 'min_mass', 'calls', 'iterations', 'warm_calls',
                                     'warm_iterations', 'seed', 'xi', 'mu_factor', 'softcut', 'collcut', 'nts')]
    args.append(config['angular_map'])
    (folder / 'configuration.json').write_text(json.dumps(config, indent=2) + '\n')
    (folder / 'build-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    start = time.monotonic()
    with (folder / 'run.log').open('w') as log:
        subprocess.run(args, cwd=folder, stdout=log, stderr=subprocess.STDOUT, check=True)
    # Fortran STOP can return zero and may leave an empty output file.
    result = json.loads((folder / 'result.json').read_text())
    if result.get('schema') != 'leptonff-mcmule-fixed-order-v1' or len(result['observables']) != 14:
        raise ValueError('Incomplete native result')
    if result['bad_phase_space_points'] != 0:
        raise ValueError('Bad phase-space points recorded')
    verify(build, manifest)
    result['wall_seconds'] = time.monotonic() - start
    result['build_manifest_sha256'] = sha(folder / 'build-manifest.json')
    (folder / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(f"{folder.name}: {result['wall_seconds']:.1f}s", flush=True)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, default=Path('build-mcmule'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--pieces', nargs='+', choices=PIECES, required=True)
    p.add_argument('--mode', choices=('shape', 'inclusive'), default='shape')
    p.add_argument('--Q', type=float, default=91.1876)
    p.add_argument('--min-mass', type=float, default=0)
    p.add_argument('--calls', type=int, default=100000)
    p.add_argument('--iterations', type=int, default=8)
    p.add_argument('--warm-calls', type=int, default=50000)
    p.add_argument('--warm-iterations', type=int, default=6)
    p.add_argument('--seed', type=int, default=831001)
    p.add_argument('--seeds', type=int, default=1)
    p.add_argument('--jobs', type=int, default=4)
    p.add_argument('--xi', type=float, default=.2)
    p.add_argument('--mu-factor', type=float, default=1)
    p.add_argument('--softcut', type=float, default=1e-10)
    p.add_argument('--collcut', type=float, default=1e-14)
    p.add_argument('--nts', type=float, default=1e-3)
    p.add_argument('--angular-map', choices=('native', 'fsr', 'fsr-direct'), default='fsr-direct')
    a = p.parse_args()
    if a.mode == 'shape' and set(a.pieces) & {'ee2mm0', 'ee2mmFMM'}:
        p.error('No-real-photon shape is zero; run these pieces in inclusive mode')
    configs = []
    for i, piece in enumerate(a.pieces):
        for j in range(a.seeds):
            seed = a.seed + 1009 * (i * a.seeds + j)
            config = {k: v for k, v in vars(a).items() if k not in ('pieces', 'jobs', 'seeds')}
            config.update(build=str(a.build.resolve()), piece=piece, seed=seed,
                          output=str(a.output.resolve() / f'{piece}-{a.mode}-{seed}'))
            configs.append(config)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(execute, configs))


if __name__ == '__main__':
    main()
