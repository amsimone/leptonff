#!/usr/bin/env python3
"""Optional native/mapped phase-space volume and kinematics cross-check."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build', type=Path, default=Path('build-mcmule'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--calls', type=int, default=1000000)
    p.add_argument('--Q', type=float, default=3)
    a = p.parse_args(); build = a.build.resolve(); here = Path(__file__).resolve().parent
    manifest = json.loads((build / 'manifest.json').read_text())
    command = manifest['command'][:]
    index = next(i for i, x in enumerate(command) if x.endswith('/reference.f90'))
    command[index] = str(here / 'check_phase_space.f90')
    command[-1] = str(build / 'check_phase_space')
    subprocess.run(command, check=True)
    rows = json.loads(subprocess.check_output([command[-1], str(a.calls), str(a.Q)], text=True))
    checks = []
    for row in rows:
        checks.append(dict(name=f"map {row['map_index']} kinematics",
                           passed=row['mass_shell_residual'] < 1e-10 and row['momentum_residual'] < 1e-10))
    for i, j in [(0, 1), (2, 3), (2, 4), (3, 4)]:
        error = math.hypot(rows[i]['mc_error'], rows[j]['mc_error'])
        pull = (rows[i]['volume']-rows[j]['volume'])/error
        checks.append(dict(name=f'volume maps {i+1} vs {j+1}', pull=pull, passed=abs(pull)<5))
    report = dict(Q=a.Q, calls_per_map=a.calls, rows=rows, checks=checks,
                  convention='Uncut Lorentz-invariant phase space; FKS photon energy factors restored')
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(checks, indent=2))
    if not all(c['passed'] for c in checks):
        sys.exit(1)


if __name__ == '__main__':
    main()
