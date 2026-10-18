#!/usr/bin/env python3
"""Validate LL/NLL routing through the compiled W+/W-/Z differential callback."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import re


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reference', type=Path, help='Same callback driver linked to the frozen pre-NLL adapter/core')
    p.add_argument('--reference-tolerance', type=float, default=0., help='Absolute tolerance for a refactor regression; zero requires bitwise equality')
    a = p.parse_args()
    checks, samples, logs = [], {}, {}
    environment = dict(os.environ, OMP_NUM_THREADS='1', LEPTONFF_PHOTONIC_LL='1',
                       LEPTONFF_ALPHA_SCHEME='msbar', LEPTONFF_FLAVOURS='one-family',
                       LEPTONFF_COUPLED_HARD='0', LEPTONFF_FRACTION_CELLS='0',
                       LEPTONFF_DRESS_HARD='1', LEPTONFF_HARD_SAMPLES='4',
                       LEPTONFF_SPECTRA='1', LEPTONFF_GRID_BINS='256', LEPTONFF_LOG_SPAN='40',
                       LEPTONFF_ALPHA=str(1/137.035999084), LEPTONFF_HARD_SOFT='1e-5',
                       LEPTONFF_HARD_COLLINEAR='1e-5', LEPTONFF_PHOTOS_ORDER='exponentiated',
                       LEPTONFF_PHOTOS_CUTOFF='1e-5', LEPTONFF_PHOTOS_LINEAR_SCALE='1',
                       LEPTONFF_PT_MIN='25', LEPTONFF_ETA_MAX='2.4', LEPTONFF_NU_PT_MIN='0',
                       LEPTONFF_MT_MIN='0', LEPTONFF_MASS_MIN='60', LEPTONFF_MASS_MAX='120',
                       LEPTONFF_THRESHOLD_FACTOR='1')
    for process in [1, 6, 31]:
        rows = []
        for nll in [False, True]:
            env = dict(environment, LEPTONFF_PHOTONIC_NLL=str(int(nll)),
                       LEPTONFF_MU_FRAG=str(91.1876 if process == 31 else 80.385))
            run = subprocess.run([str(a.executable.resolve()), str(process)], env=env,
                                 capture_output=True, text=True, check=True)
            checks.append(dict(label=f'{process} NLL={nll} applied mode',
                               passed=('leptonff photonic NLL:' in run.stderr) == nll))
            logs[f'{process}_nll_{nll}'] = run.stderr
            if nll:
                applied = re.search(r'source_N1=([^;]+); raw_inclusive_factor=([^;]+); kernel_number_residual=([^;]+); evolved_number_residual=([^;]+); wrap_bound=([^;]+); tail_bound=([^\n]+)', run.stderr)
                values = list(map(float, applied.groups())) if applied else []
                checks.append(dict(label=f'{process} explicit normalisation and numerical diagnostics',
                                   passed=len(values) == 6 and abs(values[0]+.719311605774509) < 1e-12
                                   and 0 < values[1] <= 1 and all(abs(x) < 1e-8 for x in values[2:])))
            rows.append(json.loads(run.stdout.split('LEPTONFF_JSON=', 1)[1]))
        samples[str(process)] = dict(ll=rows[0], nll=rows[1])
        if a.reference:
            env = dict(environment, LEPTONFF_PHOTONIC_NLL='0',
                       LEPTONFF_MU_FRAG=str(91.1876 if process == 31 else 80.385))
            reference = subprocess.run([str(a.reference.resolve()), str(process)], env=env,
                                       capture_output=True, text=True, check=True)
            old = json.loads(reference.stdout.split('LEPTONFF_JSON=', 1)[1])
            maximum=max(abs(x-y) for left,right in zip(old,rows[0]) for field in ['rates','spectra']
                        for x,y in zip(left[field],right[field]))
            checks.append(dict(label=f'{process} all LL rates/spectra agree with frozen adapter',
                               passed=maximum <= a.reference_tolerance,bitwise=old == rows[0],
                               maximum_absolute_difference=maximum,tolerance=a.reference_tolerance))
        changed = False
        for point, (ll, nll) in enumerate(zip(*rows)):
            unchanged_rates = [0, 1, 2, 3, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 20, 21]
            checks.append(dict(label=f'{process} point {point} common Born/NLO/hard/PHOTOS rates',
                               passed=all(ll['rates'][k] == nll['rates'][k] for k in unchanged_rates)))
            checks.append(dict(label=f'{process} point {point} spectrum Born/NLO/PHOTOS unchanged',
                               passed=all(ll['spectra'][8*k+j] == nll['spectra'][8*k+j]
                                          for k in range(17) for j in [0, 1, 5, 6])))
            changed |= abs(ll['rates'][23]-nll['rates'][23]) > 1e-10
        checks.append(dict(label=f'{process} NLL changes matched result', passed=changed))
    report = dict(scope='four phase-space points per W charge and Z, actual callback; no host integration or precision claim',
                  executable_sha256=hashlib.sha256(a.executable.read_bytes()).hexdigest(),
                  reference_executable_sha256=hashlib.sha256(a.reference.read_bytes()).hexdigest() if a.reference else None,
                  environment=environment, logs=logs,
                  checks=checks, samples=samples, passed=all(c['passed'] for c in checks))
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(checks=len(checks), failed=[c for c in checks if not c['passed']]), indent=2))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
