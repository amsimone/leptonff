#!/usr/bin/env python3
"""Independent fixed-quadrature check of the conditioned photonic NLL moments."""
import argparse
import json
import math
from pathlib import Path
import subprocess

from check_fiducial_reference import gauss_rule
from check_nlo_kernel_reference import dilog
from check_reference import source_digest


def moments(N, count=128):
    # Independently implement the published regular density; subtract the
    # endpoint analytically through (z^(N-1)-1), not a fitted delta shift.
    kernel = boundary = 0.
    for upper_half in [False, True]:
        for x, w in gauss_rule(count):
            t = (x+1)/2
            edge = .5*t**4
            z, q = (1-edge, edge) if upper_half else (edge, 1-edge)
            l, s = (math.log1p(-q), math.log(q)) if upper_half else (math.log(z), math.log1p(-z))
            test = math.expm1((N-1)*l)
            p = (1+z*z)/q
            g = p*(2*l*s-2*l*l-2*dilog(q))+(1+z)*l*l/2+2*z*l-3*z+2
            jacobian = w*t**3
            kernel += jacobian*g*test
            boundary += jacobian*p*(-1-2*s)*test
    return kernel, boundary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable', type=Path)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    data = json.loads(subprocess.check_output([str(a.executable.resolve())], text=True))
    checks = []

    def check(label, actual, expected, tolerance):
        checks.append(dict(label=label, actual=actual, expected=expected,
                           tolerance=tolerance, passed=abs(actual-expected) <= tolerance))

    aa = data['alpha']/(2*math.pi)
    L = 2*math.log(data['mu']/data['mu0'])
    source_number = -13/4+math.pi**2/2-2*1.2020569031595942854
    check('published raw first moment', data['source_number'], source_number, 2e-14)
    for row in data['rows']:
        N = int(row['N'])
        kernel, boundary = moments(N)
        p0 = 1.5-2*sum(1/k for k in range(1, N+1))+1/(N*(N+1))
        u = math.exp(aa*L*p0)
        exact = u*(1+aa*boundary+aa*aa*L*kernel)
        key = f"N={N}, bins={row['bins']}"
        check(key+' independent kernel', row['kernel'], kernel, 3e-8)
        check(key+' continuum NLL moment', row['exact'], exact, 2e-10)
        check(key+' raw inclusive factor', row['raw_inclusive_factor'], 1+aa*aa*L*source_number, 2e-14)
        if N == 1:
            check(key+' unitary number', row['moment'], 1., 2e-11)
        if row['bins'] == 65536:
            check(key+' refined transported kernel', row['transported_kernel'], u*kernel, .005*max(1., abs(u*kernel)))
            check(key+' refined complete moment', row['moment'], exact, 2e-4)
            c2 = L*L*p0*p0/2+L*(p0*boundary+kernel)
            check(key+' refined C2', row['coefficient2'], c2, .003*max(1., abs(c2)))
    for N in [2, 3, 5]:
        rows = [r for r in data['rows'] if r['N'] == N]
        first, last = rows[0], rows[-1]
        initial_error = abs(first['moment']-first['exact'])
        final_error = abs(last['moment']-last['exact'])
        checks.append(dict(label=f'N={N} grid convergence', passed=final_error < initial_error/5,
                           coarse_error=initial_error, fine_error=final_error))
    report = dict(scope='fixed-alpha unitary-conditioned photonic NLL, independent continuum moments and grid refinement',
                  source_sha256=source_digest(Path(__file__).resolve().parents[1]),
                  data=data, checks=checks, passed=all(c['passed'] for c in checks))
    if a.output:
        a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    failed = [c for c in checks if not c['passed']]
    print(json.dumps(dict(checks=len(checks), failed=failed), indent=2))
    if failed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
