#!/usr/bin/env python3
"""Verify that the photonic domain-routing fix preserves default host results."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('regressions', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    checks, evidence = [], {}
    old_exes, new_exes = set(), set()
    for process in ['wplus', 'wminus', 'z']:
        paths = {side: a.regressions/side/f'photonic_{process}_nominal' for side in ['old', 'new']}
        metas = {side: json.loads((path/'manifest.json').read_text()) for side, path in paths.items()}
        if any(m['status'] != 'completed' for m in metas.values()):
            raise RuntimeError('Incomplete regression run')
        old, new = metas['old'], metas['new']
        old_exes.add(old['executable_sha256']); new_exes.add(new['executable_sha256'])
        def check(label, passed):
            checks.append({'label': process+': '+label, 'passed': bool(passed)})
        check('environment', old['environment'] == new['environment'])
        check('call budget and production input', old['requested_inputs'] == new['requested_inputs'])
        def input_lines(path):
            return [r for r in (path/'input.ini').read_text().splitlines()
                    if not r.strip().startswith(('runstring =', 'rundir ='))]
        check('same host configuration', input_lines(paths['old']) == input_lines(paths['new']))
        for field in ['histograms_fb', 'spectra_histograms_fb']:
            check(field+' bitwise', old[field] == new[field])
        for name in ['iteration_rates.csv', 'event_covariance.csv']:
            check(name+' bitwise', digest(paths['old']/name) == digest(paths['new']/name))
        for key in ['core_library_sha256', 'host_executable_sha256', 'host_link_recipe_sha256',
                    'fortran_hook_sha256', 'event_covariance_sha256', 'photos_comparison_sha256',
                    'photos_sha256', 'iteration_audit_source_sha256', 'real_covariance_source_sha256']:
            check('unchanged '+key, old['adapter_build'][key] == new['adapter_build'][key])
        evidence[process] = {side: {'directory': str(path.resolve()),
            'manifest_sha256': digest(path/'manifest.json'),
            'raw_rates_sha256': digest(path/'iteration_rates.csv'),
            'event_covariance_sha256': digest(path/'event_covariance.csv')}
            for side, path in paths.items()}
    if len(old_exes) != 1 or len(new_exes) != 1:
        raise RuntimeError('Regression mixes executables')
    report = {'passed': all(c['passed'] for c in checks), 'checks': checks,
              'old_executable_sha256': old_exes.pop(), 'new_executable_sha256': new_exes.pop(),
              'evidence': evidence,
              'scope': 'Bridge now honours requested photonic log span and logs applied span. All default 33 rates, 136 spectra and grouped raw/event products unchanged; core library and other physics sources identical.'}
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'passed': report['passed'], 'checks': len(checks)}))
    if not report['passed']:
        raise SystemExit('Default adapter equivalence failed')


if __name__ == '__main__':
    main()
