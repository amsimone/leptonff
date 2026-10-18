#!/usr/bin/env python3
"""Run pointwise DYTurbo W NLO+NLL decay coupling and export reusable QCD records."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
from run import Mode, sha256, steering


def decay_steering(*, process='wp', calls=10000, seed=123456):
    if process not in ('wp', 'wm') or calls < 100 or calls % 10:
        raise ValueError('Requires W+/W- and a call budget divisible by ten')
    return steering(Mode('nlo-nll-vegas', 1, False, 'vegas'), process=process,
                    mass=(60, 100), calls=calls, seed=seed) + (
        'finord = false\nintDimRes = 4\nintDimCT = 6\nintDimVJ = 7\n'
        'mellin1d = false\nptl_bins = [ 0 1000000 ]\netal_bins = [ -100 100 ]\n')


def ratio(mean, covariance, numerator, denominator):
    if mean[denominator] <= 0:
        return {'value': None, 'error': None}
    r = mean[numerator] / mean[denominator]
    v = (covariance[5*numerator+numerator] + r*r*covariance[5*denominator+denominator]
         - 2*r*covariance[5*numerator+denominator]) / mean[denominator]**2
    if v < -1e-10:
        raise ValueError('Negative ratio variance')
    return {'value': r, 'error': math.sqrt(max(0, v))}


def assess(data, target=.001):
    if data['schema'] != 'leptonff-dyturbo-decay-1' or len(data['bins']) != 49:
        raise ValueError('Unexpected differential decay result')
    out = []
    for b in data['bins']:
        m, c = b['mean'], b['covariance']
        if len(m) != 5 or len(c) != 25 or not all(map(math.isfinite, m+c)):
            raise ValueError('Invalid differential result or covariance')
        if any(not math.isclose(c[5*a+d], c[5*d+a], rel_tol=1e-9, abs_tol=1e-6)
               for a in range(5) for d in range(5)):
            raise ValueError('Asymmetric covariance')
        out.append({'lower': b['lower'], 'upper': b['upper'],
                    'order_alpha': ratio(m, c, 1, 0),
                    'photonic_nll': ratio(m, c, 2, 0), 'photos': ratio(m, c, 4, 3)})
    errors = [b[k]['error'] for b in out[1:] for k in ('photonic_nll', 'photos')]
    return {'target_absolute_ratio_error': target,
            'precision_met': all(e is not None and e <= target for e in errors),
            'ratios': out,
            'note': 'Integration errors only; grid, geometry, cutoffs, QCD controls and independent seeds require separate checks.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', type=Path, default=Path('build-dyturbo-repro/leptonff_dyturbo_decay'))
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--pdf-data', required=True, type=Path)
    p.add_argument('--process', choices=('wp', 'wm'), default='wp')
    p.add_argument('--calls', type=int, default=10000)
    p.add_argument('--seed', type=int, default=123456)
    p.add_argument('--grid', type=int, default=8192)
    p.add_argument('--hard-samples', type=int, default=8)
    p.add_argument('--photos-samples', type=int, default=32)
    p.add_argument('--alpha', type=float, default=1/137.035999084)
    a = p.parse_args()
    card = decay_steering(process=a.process, calls=a.calls, seed=a.seed)
    exe, out = a.executable.resolve(), a.output.resolve()
    backends = (exe.parent/'decay-backend-paths.txt').read_text().splitlines()
    provenance = {'client_sha256': sha256(exe),
                  'backend_sha256': [sha256(path) for path in backends],
                  'input_sha256': hashlib.sha256(card.encode()).hexdigest()}
    grid = a.pdf_data.resolve()/'NNPDF31_nnlo_as_0118'
    provenance['pdf_sha256'] = [sha256(grid/name) for name in
                              ('NNPDF31_nnlo_as_0118.info', 'NNPDF31_nnlo_as_0118_0000.dat')]
    identity = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True, exist_ok=False)
    (out/'input.in').write_text(card)
    env = dict(os.environ, LHAPDF_DATA_PATH=str(a.pdf_data.resolve()), OMP_NUM_THREADS='1',
               LEPTONFF_DYTURBO_RESULT=str(out/'native.json'), LEPTONFF_RECORDS=str(out/'production.lff'),
               LEPTONFF_PRODUCTION_ID=identity, LEPTONFF_CALLS=str(a.calls), LEPTONFF_ALPHA=str(a.alpha),
               LEPTONFF_GRID_BINS=str(a.grid), LEPTONFF_HARD_SAMPLES=str(a.hard_samples),
               LEPTONFF_PHOTOS_SAMPLES=str(a.photos_samples))
    (out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    with (out/'backend.log').open('w') as log:
        subprocess.run([str(exe), str(out/'input.in')], cwd=out, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    data = json.loads((out/'native.json').read_text())
    assessment = assess(data)
    (out/'precision.json').write_text(json.dumps(assessment, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'output': str(out), 'precision_met': assessment['precision_met'],
                      'total_ratios': assessment['ratios'][0]}, indent=2))


if __name__ == '__main__':
    main()
