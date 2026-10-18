#!/usr/bin/env python3
"""Check native-client closure against DYTurbo's original executable."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess

from run import sha256


def compare(directory, executable, *, pdf_data=None, timeout=1800):
    directory, executable = Path(directory).resolve(), Path(executable).resolve()
    result = json.loads((directory / 'result.json').read_text())
    reference = directory / 'reference'
    reference.mkdir(exist_ok=False)
    card = (directory / 'input.in').read_text()
    if card.count('silent = true\n') != 1:
        raise ValueError('Unrecognised native card')
    (reference / 'input.in').write_text(card.replace('silent = true\n', 'silent = false\n'))
    env = dict(os.environ, OMP_NUM_THREADS='1', CUBACORES='0')
    data = pdf_data or result['provenance']['pdf_data_path']
    if data:
        env['LHAPDF_DATA_PATH'] = str(data)
    with (reference / 'backend.log').open('w') as log:
        subprocess.run([str(executable), 'input.in'], cwd=reference, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=timeout)
    rows = [list(map(float, line.split())) for line in (reference / 'results.txt').read_text().splitlines()
            if line.strip() and not line.startswith('#')]
    if len(rows) != len(result['bins']) + 1:
        raise ValueError('Unexpected original DYTurbo result table')
    for b, row in zip(result['bins'], rows):
        if len(row) != 2 or not all(math.isfinite(x) for x in row):
            raise ValueError('Invalid reference output')
        for key, expected in zip(('value', 'error'), row):
            if not math.isclose(b[key], expected, rel_tol=1e-11, abs_tol=1e-7):
                raise ValueError(f"Original-executable mismatch: {key} {b[key]} vs {expected}")
    report = dict(status='passed', reference_executable_sha256=sha256(executable),
                  scope='Same-seed QCD-only native library client vs original executable; no QED closure claim',
                  bins=len(result['bins']), relative_tolerance=1e-11, absolute_tolerance_fb=1e-7)
    (reference / 'closure.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True, help='Matrix output directory')
    p.add_argument('--executable', type=Path, required=True, help='Original DYTurbo executable')
    p.add_argument('--pdf-data', type=Path)
    p.add_argument('--timeout', type=float, default=1800)
    a = p.parse_args()
    reports = {}
    for path in sorted(a.results.glob('*/result.json')):
        print('Checking ' + path.parent.name, flush=True)
        reports[path.parent.name] = compare(path.parent, a.executable, pdf_data=a.pdf_data, timeout=a.timeout)
        print('  original-executable closure passed', flush=True)
    if not reports:
        p.error('No native results found')
    (a.results / 'reference_closure.json').write_text(json.dumps(reports, indent=2) + '\n')


if __name__ == '__main__':
    main()
