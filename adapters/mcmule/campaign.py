#!/usr/bin/env python3
"""Optional fixed-order logarithm scan or FKS/technical-control campaign."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from run import execute


def configurations(build, output, kind):
    points = ([dict(Q=q, xi=.2, softcut=1e-10, collcut=1e-14, nts=1e-3)
               for q in (10., 30., 60., 300., 1000.)] if kind == 'log-scan' else
              [dict(Q=91.1876, xi=.05, softcut=1e-10, collcut=1e-14, nts=1e-3),
               dict(Q=91.1876, xi=.5, softcut=1e-10, collcut=1e-14, nts=1e-3),
               dict(Q=91.1876, xi=.2, softcut=1e-11, collcut=1e-15, nts=1e-4)])
    configs = []
    stream = 1001001 if kind == 'log-scan' else 2001001
    for point in points:
        label = (f"Q{point['Q']:g}" if kind == 'log-scan' else
                 f"xi{point['xi']:g}-soft{point['softcut']:g}")
        for mode, piece in [('inclusive', 'ee2mmFMM'), ('inclusive', 'ee2mmRMM'),
                            ('shape', 'ee2mmRMM'), ('shape', 'ee2mmRRMMMM'), ('shape', 'ee2mmRFMMMM')]:
            # Two independent streams per expensive/singular shape piece.
            for replica in range(2 if mode == 'shape' else 1):
                stream += 1009
                rv = piece == 'ee2mmRFMMMM'
                config = dict(build=str(build.resolve()), mode=mode, piece=piece, min_mass=0,
                              calls=25000 if rv else 1000000, iterations=8,
                              warm_calls=5000 if rv else 100000, warm_iterations=6,
                              seed=stream, mu_factor=1, angular_map='fsr-direct', **point)
                config['output'] = str(output.resolve() / label / f'{piece}-{mode}-{stream}')
                configs.append(config)
    return configs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('kind', choices=('log-scan', 'controls'))
    p.add_argument('--build', type=Path, default=Path('build-mcmule'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--jobs', type=int, default=6)
    a = p.parse_args(); configs = configurations(a.build, a.output, a.kind)
    a.output.mkdir(parents=True, exist_ok=False)
    (a.output / 'campaign.json').write_text(json.dumps(configs, indent=2) + '\n')
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(execute, configs))


if __name__ == '__main__':
    main()
