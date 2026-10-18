#!/usr/bin/env python3
"""Run native DYTurbo QCD production with explicit accuracy and integration settings."""
import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


@dataclass(frozen=True)
class Mode:
    name: str
    order: int
    fixed: bool
    integration: str


MODES = {m.name: m for m in (
    Mode('nlo-quad', 1, True, 'quadrature'),
    Mode('nlo-vegas', 1, True, 'vegas'),
    Mode('nnlo-quad', 2, True, 'quadrature'),
    Mode('nnlo-vegas', 2, True, 'vegas'),
    Mode('nlo-nll', 1, False, 'quadrature'),
    Mode('nnlo-nnll', 2, False, 'quadrature'),
)}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def steering(mode, *, process='z', pdf='NNPDF31_nnlo_as_0118', member=0,
             mass=(66., 116.), qt=(0., 6500.), rapidity=(-10., 10.),
             calls=100000, iterations=5, accuracy=.001, seed=123456,
             fiducial=False, hybrid=False, mass_bin_width=None):
    """Return a complete override card; the backend supplies its versioned defaults."""
    if process not in ('wp', 'wm', 'z'):
        raise ValueError('Unsupported Drell-Yan process')
    if mode.order not in (1, 2) or mode.integration not in ('quadrature', 'vegas'):
        raise ValueError('Unsupported production accuracy or integration method')
    if any(len(b) != 2 or not all(math.isfinite(x) for x in b) or b[0] >= b[1]
           for b in (mass, qt, rapidity)) or mass[0] <= 0 or qt[0] < 0:
        raise ValueError('Invalid production bounds')
    if any(type(x) is not int or x < 0 for x in (member, seed, iterations, calls)) or calls < 100:
        raise ValueError('Invalid sampling settings')
    if not math.isfinite(accuracy) or not 0 < accuracy < 1:
        raise ValueError('Invalid requested relative accuracy')
    if not pdf or any(c.isspace() or c in '#=[]' for c in pdf):
        raise ValueError('Invalid PDF name')
    quad = mode.integration == 'quadrature'
    if fiducial and quad and mode.order == 2 and not hybrid:
        raise ValueError('DYTurbo 1.4.2 cannot integrate NNLO V+jet with lepton cuts by quadrature; '
                         'use VEGAS or explicitly request --hybrid')
    settings = {
        'nproc': {'wp': 1, 'wm': 2, 'z': 3}[process],
        'sroot': 13000, 'LHAPDFset': pdf, 'LHAPDFmember': member,
        'fixedorder_only': mode.fixed, 'order': mode.order, 'primed': True,
        # Keep the hard/CT/V+jet matching at the requested QCD order. In
        # 1.4.2 unprimed defaults lower these orders. Record the convention.
        'order_sudak': mode.order, 'order_hcoef': mode.order,
        'order_evol': mode.order, 'order_expc': mode.order,
        'order_ct': mode.order, 'order_vjet': mode.order, 'orderlha': False,
        'qed': False, 'makecuts': fiducial, 'makehistos': False,
        'BORNquad': quad, 'CTquad': quad,
        'VJquad': quad and not (fiducial and mode.order == 2), 'FPCquad': quad,
        'cores': 0, 'rseed': seed, 'silent': True, 'redirect': False,
        'niterBORN': iterations, 'niterCT': iterations, 'niterVJ': iterations,
        'vegasncallsBORN': calls, 'vegasncallsCT': calls,
        'vegasncallsVJLO': calls, 'vegasncallsVJREAL': calls,
        'vegasncallsVJVIRT': calls, 'relaccuracy': accuracy,
        'pcubature': True, 'ptbinwidth': False, 'ybinwidth': False,
        'mbinwidth': False, 'dsdqt2': False, 'edsdp3': False,
        'force_binsampling': True,
        'm_bins': mass, 'qt_bins': qt, 'y_bins': rapidity,
    }
    if mass_bin_width is not None:
        if not math.isfinite(mass_bin_width) or mass_bin_width <= 0:
            raise ValueError('Invalid production mass bin width')
        count = round((mass[1]-mass[0])/mass_bin_width)
        if not 1 <= count <= 10000 or not math.isclose(count*mass_bin_width,mass[1]-mass[0],rel_tol=1e-12):
            raise ValueError('Mass bin width must partition the requested mass interval')
        settings['m_bins'] = [mass[0]+i*mass_bin_width for i in range(count+1)]
    def value(x):
        if isinstance(x, bool):
            return str(x).lower()
        if isinstance(x, (tuple, list)):
            return '[ ' + ' '.join(format(y, '.17g') for y in x) + ' ]'
        return str(x)
    return ''.join(f'{k} = {value(v)}\n' for k, v in settings.items())


def validate_result(result, mode, *, fiducial=False, hybrid=False):
    if result.get('schema') != 'leptonff-dyturbo-native-1' or result.get('units') != 'fb':
        raise ValueError('Unexpected native result schema')
    if result['qcd_order'] != mode.order or result['fixed_order'] != mode.fixed:
        raise ValueError('Backend did not resolve the requested QCD prediction')
    if result['fiducial'] != fiducial or not result['bins']:
        raise ValueError('Backend resolved different cuts or no integration bins')
    orders = result['resolved_orders']
    if any(orders[k] != mode.order for k in ('hard', 'counterterm', 'vjet', 'sudakov')):
        raise ValueError('Backend lowered the perturbative matching order')
    for b in result['bins']:
        terms = b['terms']
        if not terms or abs(sum(t['value'] for t in terms) - b['value']) > 1e-10 * max(1., abs(b['value'])):
            raise ValueError('Native term sum is inconsistent')
        for item in [b] + terms:
            if not all(math.isfinite(item[k]) for k in ('value', 'error')) or item['error'] < 0:
                raise ValueError('Invalid native cross section or uncertainty')
        if mode.integration == 'vegas' and any(t['method'] != 'vegas' for t in terms):
            raise ValueError('Requested VEGAS, but the backend used quadrature')
        if mode.integration == 'quadrature' and not hybrid and any(t['method'] != 'quadrature' for t in terms):
            raise ValueError('Requested quadrature, but the backend used VEGAS')
        # Error combination is DYTurbo's term-wise prescription, not a
        # covariance inferred by treating quadrature nodes as random events.
        if not math.isclose(math.hypot(*(t['error'] for t in terms)), b['error'], rel_tol=1e-12):
            raise ValueError('Native error combination is inconsistent')


def run(executable, directory, mode, *, pdf_data=None, timeout=1800, **settings):
    executable, directory = Path(executable).resolve(), Path(directory).resolve()
    provenance = {'executable_sha256': sha256(executable)}
    backend_paths = executable.parent / 'backend-paths.txt'
    if backend_paths.is_file():
        library, defaults = backend_paths.read_text().splitlines()
        provenance.update(backend_library_sha256=sha256(library), default_card_sha256=sha256(defaults))
    if pdf_data:
        grid = Path(pdf_data).resolve() / settings.get('pdf', 'NNPDF31_nnlo_as_0118')
        info = grid / (grid.name + '.info')
        member = grid / f"{grid.name}_{settings.get('member', 0):04d}.dat"
        provenance.update(pdf_info_sha256=sha256(info), pdf_member_sha256=sha256(member))
    directory.mkdir(parents=True, exist_ok=False)
    card = directory / 'input.in'
    card.write_text(steering(mode, **settings))
    env = dict(os.environ, LEPTONFF_DYTURBO_RESULT=str(directory / 'native.json'),
               OMP_NUM_THREADS='1', CUBACORES='0')
    if pdf_data:
        env['LHAPDF_DATA_PATH'] = str(Path(pdf_data).resolve())
    start = time.monotonic()
    with (directory / 'backend.log').open('w') as log:
        subprocess.run([str(executable), str(card)], cwd=directory, env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=timeout)
    result = json.loads((directory / 'native.json').read_text())
    validate_result(result, mode, fiducial=settings.get('fiducial', False), hybrid=settings.get('hybrid', False))
    result['requested_mode'] = mode.name
    provenance.update(input_sha256=sha256(card), pdf_data_path=env.get('LHAPDF_DATA_PATH'),
                      wall_seconds=time.monotonic() - start)
    result['provenance'] = provenance
    result['accuracy_note'] = ('Pilot or finite-budget integration. Reported errors are backend estimates; '
                              'quadrature refinement, independent VEGAS seeds and qt-cutoff scans are required '
                              'before a phenomenology accuracy claim.')
    (directory / 'result.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--executable', type=Path, default=Path('build-dyturbo-interface/leptonff_dyturbo'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mode', choices=list(MODES), action='append')
    p.add_argument('--process', choices=['z', 'wp', 'wm'], default='z')
    p.add_argument('--pdf', default='NNPDF31_nnlo_as_0118')
    p.add_argument('--member', type=int, default=0)
    p.add_argument('--pdf-data', type=Path)
    p.add_argument('--mass', type=float, nargs=2, default=[66., 116.])
    p.add_argument('--mass-bin-width', type=float, help='Partition parent-mass range for reusable mass migration')
    p.add_argument('--qt', type=float, nargs=2, default=[0., 6500.])
    p.add_argument('--rapidity', type=float, nargs=2, default=[-10., 10.])
    p.add_argument('--calls', type=int, default=100000)
    p.add_argument('--iterations', type=int, default=5)
    p.add_argument('--accuracy', type=float, default=.001)
    p.add_argument('--seed', type=int, default=123456)
    p.add_argument('--timeout', type=float, default=1800)
    p.add_argument('--fiducial', action='store_true', help='Native Born-lepton cuts, without QED decay')
    p.add_argument('--hybrid', action='store_true', help='Explicitly allow NNLO fiducial V+jet VEGAS in a quadrature run')
    a = p.parse_args()
    modes = [MODES[x] for x in (a.mode or list(MODES))]
    settings = {k: getattr(a, k) for k in ('process', 'pdf', 'member', 'mass', 'qt', 'rapidity',
                                          'calls', 'iterations', 'accuracy', 'seed', 'fiducial', 'hybrid', 'mass_bin_width')}
    # Validate the entire matrix before starting any production run.
    for mode in modes:
        steering(mode, **settings)
    a.output.mkdir(parents=True, exist_ok=False)
    results = []
    for mode in modes:
        print(f'Running {a.process} {mode.name} ...', flush=True)
        r = run(a.executable, a.output / mode.name, mode, pdf_data=a.pdf_data, timeout=a.timeout, **settings)
        results.append(r)
        b = r['bins'][0]
        print(f"  {b['value'] / 1000:.8g} +/- {b['error'] / 1000:.3g} pb; {r['seconds']:.2f} s", flush=True)
    (a.output / 'matrix.json').write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
