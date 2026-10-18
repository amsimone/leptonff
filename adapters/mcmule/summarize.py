#!/usr/bin/env python3
"""Report independent fixed-order logarithm fits and subtraction controls."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import chi2
from analyse import analyse


def fit(L, values, errors, powers, power_control=None):
    L = np.asarray(L); values = np.asarray(values); errors = np.asarray(errors)
    reference = 2*math.log(91.1876/.1056583755)
    X = np.column_stack([(L-reference)**p for p in powers])
    if power_control is not None:
        X = np.column_stack([X, power_control])
    covariance = np.linalg.inv((X/errors[:, None]).T @ (X/errors[:, None]))
    coeff = covariance @ (X.T @ (values/errors**2))
    pulls = (values-X@coeff)/errors
    statistic = float(pulls @ pulls); dof = len(values)-len(coeff)
    return dict(powers=powers, reference_L=reference, coefficients=coeff.tolist(),
                errors=np.sqrt(np.diag(covariance)).tolist(), covariance=covariance.tolist(),
                chi2=statistic, dof=dof, p_value=float(chi2.sf(statistic, dof)) if dof else None,
                pulls=pulls.tolist(), power_control=power_control is not None)


def load_group(folder):
    paths = sorted(folder.rglob('result.json'))
    return analyse([json.loads(p.read_text()) for p in paths])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path('results/mcmule-fixed-order'))
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); root = a.root; a.output.mkdir(parents=True, exist_ok=True)
    central = json.loads((root/'central-report.json').read_text())
    scans = [load_group(folder) for folder in sorted((root/'log-scan').glob('Q*')) if folder.is_dir()]
    scans.append(central); scans.sort(key=lambda r: r['settings']['Q'])
    controls = [load_group(folder) for folder in sorted((root/'controls').iterdir()) if folder.is_dir()]
    for group in scans:
        if any(group['settings'][k] != central['settings'][k]
               for k in central['settings'] if k not in ('Q', 'born_alpha_one_MeV_minus2')):
            raise ValueError('Log-scan changes settings beyond Q')
    for group in controls:
        if any(group['settings'][k] != central['settings'][k]
               for k in central['settings'] if k not in ('xi', 'softcut', 'collcut', 'nts_switch')):
            raise ValueError('Technical control changes physical settings')
    per_moment = []
    for i in (1, 2, 3, 4):
        points = [dict(Q=r['settings']['Q'], L=2*math.log(r['settings']['Q']/r['settings']['mass']),
                       delta_ppm=r['rows'][i]['nnlo_delta_ppm'], error_ppm=r['rows'][i]['nnlo_mc_error_ppm'],
                       nlo_delta_ppm=r['rows'][i]['nlo_delta_ppm'], nlo_error_ppm=r['rows'][i]['nlo_mc_error_ppm'])
                  for r in scans]
        def fitted(selected, powers, power=False):
            return fit([x['L'] for x in selected], [x['delta_ppm'] for x in selected],
                       [x['error_ppm'] for x in selected], powers,
                       [(10/x['Q'])**2 for x in selected] if power else None)
        per_moment.append(dict(N=central['rows'][i]['value'], points=points,
                               constant=fitted(points, [0]), linear=fitted(points, [0, 1]),
                               quadratic=fitted(points, [0, 1, 2]),
                               high_Q_linear=fitted([x for x in points if x['Q']>=30], [0, 1]),
                               linear_power_control=fitted(points, [0, 1], True)))
    comparisons = []
    for control in controls:
        rows = []
        for c, n in zip(control['rows'], central['rows']):
            error = math.hypot(c['nnlo_mc_error_ppm'], n['nnlo_mc_error_ppm'])
            rows.append(dict(observable=c['observable'], value=c['value'],
                             difference_ppm=c['nnlo_delta_ppm']-n['nnlo_delta_ppm'],
                             mc_error_ppm=error,
                             pull=(c['nnlo_delta_ppm']-n['nnlo_delta_ppm'])/error if error else 0))
        comparisons.append(dict(settings=control['settings'], rows=rows))
    report = dict(schema='leptonff-mcmule-summary-v1', central=central, scan=scans,
                  log_fits=per_moment, controls=comparisons,
                  convention='Absolute normalized-moment/cumulative ppm; alpha^2 coefficient difference only',
                  caveats=['Fits are individual-observable fits, not a global correlated chi-square.',
                           'The inverse-Q-squared fit is a power-sensitivity diagnostic, not a derived massive remainder.',
                           'Errors are independent-piece/seed MC estimates. Fits do not rescale covariance by chi-square.',
                           'No production kernel is fitted or changed.'])
    paths = [Path(item['path']) for item in central['inputs']]
    paths += sorted((root/'log-scan').rglob('result.json'))
    paths += sorted((root/'controls').rglob('result.json'))
    if len(set(x.resolve() for x in paths)) != len(paths):
        raise ValueError('Repeated raw input')
    used_seeds = set()
    for path in paths:
        seed = json.loads(path.read_text())['seed']
        if seed in used_seeds:
            raise ValueError('Scan/control stream reused; cross-point covariance unavailable')
        used_seeds.add(seed)
    report['inputs'] = [dict(path=str(x), sha256=hashlib.sha256(x.read_bytes()).hexdigest()) for x in paths]
    here = Path(__file__).resolve().parent
    sources = [here/'summarize.py', here/'analyse.py', here.parent/'validation/check_mass_radiator.py']
    report['analysis_sources'] = {str(x): hashlib.sha256(x.read_bytes()).hexdigest() for x in sources}
    (a.output/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for axis, row in zip(axes.flat, per_moment):
        points = row['points']; constant = row['constant']
        axis.errorbar([x['L'] for x in points], [x['delta_ppm'] for x in points],
                      yerr=[x['error_ppm'] for x in points], fmt='o', label='McMule − leptonff NLL expansion')
        c = constant['coefficients'][0]; e = constant['errors'][0]
        axis.axhline(c, color='C1', label='Constant fit')
        axis.axhspan(c-e, c+e, color='C1', alpha=.2)
        axis.set_title(f"N={row['N']:g}")
        axis.set_xlabel('L = log(Q²/mμ²)'); axis.set_ylabel('α² contribution difference [ppm]')
        axis.grid(alpha=.2)
    axes.flat[0].legend(fontsize=7)
    fig.suptitle('Pure vector current, bare muons, fixed α; normalized photonic FSR')
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(a.output/('fixed-order-log-scan.'+suffix), dpi=180)
    plt.close(fig)
    for row in per_moment:
        c, l = row['constant'], row['high_Q_linear']
        print(f"N={row['N']:g}: constant {c['coefficients'][0]:.3f} +/- {c['errors'][0]:.3f} ppm; "
              f"chi²/dof {c['chi2']:.2f}/{c['dof']}; high-Q slope {l['coefficients'][1]:.3f} +/- {l['errors'][1]:.3f} ppm/L")
    print('Largest control pull:', max(abs(x['pull']) for c in comparisons for x in c['rows']))


if __name__ == '__main__':
    main()
