#!/usr/bin/env python3
"""Combine fixed-order shape terms with the required inclusive NLO counterterm."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'adapters/validation'))
from check_mass_radiator import prediction, cumulative_prediction

SHAPE_PIECES = ('ee2mmRMM', 'ee2mmRFMMMM', 'ee2mmRRMMMM')
NORMALIZATION_PIECES = ('ee2mmFMM', 'ee2mmRMM')
COMMON = ('schema', 'Q', 'mass', 'electron_mass', 'support_min_mass', 'physical_alpha',
          'internal_alpha', 'born_alpha_one_MeV_minus2', 'xi', 'mu_factor', 'softcut', 'collcut',
          'current', 'fsr', 'isr', 'ifi', 'pairs', 'weak', 'coupling', 'angular_map')


def validate(records):
    if not records:
        raise ValueError('No complete results')
    base = records[0]; used = set()
    radiative_nts = {r.get('nts_switch') for r in records
                     if r['piece'] not in ('ee2mm0', 'ee2mmFMM')}
    if len(radiative_nts) != 1:
        raise ValueError('Different near-to-soft settings; keep controls separate')
    for record in records:
        if any(record[k] != base[k] for k in COMMON):
            raise ValueError('Different physics, subtraction settings or importance map')
        if record['schema'] != 'leptonff-mcmule-fixed-order-v1' or record['internal_alpha'] != 1:
            raise ValueError('Unsupported result/coupling normalization')
        if (not record['fsr'] or any(record[k] for k in ('isr', 'ifi', 'pairs', 'weak'))
                or record.get('running_alpha', False) or record['current'] != 'pure vector'):
            raise ValueError('Not the fixed-alpha pure-vector photonic FSR sector')
        if record['piece'] not in ('ee2mm0', 'ee2mmFMM') + SHAPE_PIECES:
            raise ValueError('Unexpected perturbative piece')
        if record['support_min_mass'] != 0:
            raise ValueError('Conditional support needs its own FF coefficients; not implemented here')
        if record['bad_phase_space_points'] != 0:
            raise ValueError('Bad phase-space points')
        if record.get('seed') in used:
            raise ValueError('Reused random stream; independent-piece covariance unavailable')
        used.add(record['seed'])
        if len(record['observables']) != 14:
            raise ValueError('Missing observables')
        if record['measurement'] == 'shape':
            first = record['observables'][0]
            if first['integral'] != 0 or first['mc_error'] != 0:
                raise ValueError('F-1 must vanish exactly for the normalization moment')
        elif record['measurement'] == 'inclusive':
            first = record['observables'][0]
            if any(o['integral'] != first['integral'] or o['mc_error'] != first['mc_error']
                   for o in record['observables']):
                raise ValueError('Inclusive measurement differs between observables')
        else:
            raise ValueError('Unknown measurement mode')
        for a, b in zip(record['observables'], base['observables']):
            if (a['observable'], a['value']) != (b['observable'], b['value']):
                raise ValueError('Different observable definitions')
            if not (math.isfinite(a['integral']) and math.isfinite(a['mc_error']) and a['mc_error'] >= 0):
                raise ValueError('Invalid moment or uncertainty')
        # Same release and exact source content; client hashes may differ only
        # if the caller has established and documented their equivalence.
        if record.get('build_manifest_sha256') != base.get('build_manifest_sha256'):
            raise ValueError('Different build manifests; keep controls separate')
    return base


def combine(records, piece, mode, index):
    selected = [r for r in records if r['piece'] == piece and r['measurement'] == mode]
    if not selected:
        raise ValueError('Missing piece: ' + piece + ' / ' + mode)
    # Equal independent-seed weights avoid weighting by a noisy VEGAS variance.
    values = [r['observables'][index]['integral'] for r in selected]
    mean = sum(values) / len(values)
    variance = sum(r['observables'][index]['mc_error']**2 for r in selected) / len(values)**2
    scatter = (sum((v - mean)**2 for v in values) / (len(values) * (len(values) - 1))
               if len(values) > 1 else None)
    return dict(value=mean, variance=variance, seed_scatter_variance=scatter, seeds=len(values))


def normalized_coefficients(I1, variance1, I2, variance2, c1, variance_c1):
    """All three inputs are independently integrated coefficients in alpha."""
    # Include variance of the independent product, not just its linearization.
    return dict(c1=I1, variance1=variance1, c2=I2-c1*I1,
                variance2=variance2+c1*c1*variance1+I1*I1*variance_c1+variance_c1*variance1)


def analyse(records):
    base = validate(records)
    born = base['born_alpha_one_MeV_minus2']; alpha = base['physical_alpha']
    f = combine(records, 'ee2mmFMM', 'inclusive', 0)
    r = combine(records, 'ee2mmRMM', 'inclusive', 0)
    c1 = (f['value'] + r['value']) / born
    vc1 = (f['variance'] + r['variance']) / born**2
    rows = []
    for i, observable in enumerate(base['observables']):
        r1 = combine(records, 'ee2mmRMM', 'shape', i)
        rv = combine(records, 'ee2mmRFMMMM', 'shape', i)
        rr = combine(records, 'ee2mmRRMMMM', 'shape', i)
        coeff = normalized_coefficients(r1['value']/born, r1['variance']/born**2,
                                       (rv['value']+rr['value'])/born,
                                       (rv['variance']+rr['variance'])/born**2, c1, vc1)
        args = dict(Q=base['Q'], mass=base['mass'], scale=base['Q'], alpha=alpha)
        if observable['observable'] == 'moment':
            ff = prediction(dict(args, N=observable['value']))
        else:
            ff = cumulative_prediction(dict(args, cut=observable['value']))
            ff['r1'] = (ff['order_alpha']-1) * (2*math.pi)/alpha
        expected1 = ff['r1']/(2*math.pi)
        expected2 = ff['nll2']/(2*math.pi)**2
        row = dict(observable=observable['observable'], value=observable['value'],
                   coefficient_alpha=coeff['c1'], coefficient_alpha_mc_error=math.sqrt(coeff['variance1']),
                   coefficient_alpha2=coeff['c2'], coefficient_alpha2_mc_error=math.sqrt(coeff['variance2']),
                   ff_coefficient_alpha=expected1, ff_coefficient_alpha2_nll=expected2,
                   nlo_delta_ppm=1e6*alpha*(coeff['c1']-expected1),
                   nlo_mc_error_ppm=1e6*alpha*math.sqrt(coeff['variance1']),
                   nnlo_delta_ppm=1e6*alpha**2*(coeff['c2']-expected2),
                   nnlo_mc_error_ppm=1e6*alpha**2*math.sqrt(coeff['variance2']),
                   nlo=1+alpha*coeff['c1'], nnlo=1+alpha*coeff['c1']+alpha**2*coeff['c2'],
                   pieces=dict(real=r1, real_virtual=rv, double_real=rr))
        rows.append(row)
    settings = {k: base[k] for k in COMMON}
    settings['nts_switch'] = next(r['nts_switch'] for r in records if r['piece'] == 'ee2mmRFMMMM')
    return dict(schema='leptonff-mcmule-fixed-order-comparison-v1',
                scope='Normalized pure-vector photonic FSR; not full Z NNLO electroweak',
                settings=settings,
                inclusive_nlo=dict(coefficient_alpha=c1, mc_error=math.sqrt(vc1),
                                   leading_power=3/(4*math.pi)),
                second_order='I2 - inclusive_c1 * I1; double-virtual cancels on F-1',
                uncertainty='Independent-piece VEGAS MC errors; no cross-observable covariance',
                interpretation='NLL FF omits finite alpha^2 constants and mass powers. A nonzero single-Q difference alone does not establish a logarithmic mismatch.',
                rows=rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('inputs', nargs='+', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); paths = []
    for path in a.inputs:
        paths.extend(sorted(path.rglob('result.json')) if path.is_dir() else [path])
    if len(set(paths)) != len(paths):
        raise ValueError('Repeated input path')
    records = [json.loads(path.read_text()) for path in paths]
    report = analyse(records)
    report['inputs'] = [dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()) for path in paths]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(report, indent=2) + '\n')
    for row in report['rows']:
        if row['observable'] == 'moment':
            print(f"N={row['value']:g}: NLO delta={row['nlo_delta_ppm']:.3f} +/- {row['nlo_mc_error_ppm']:.3f} ppm; "
                  f"NNLO delta={row['nnlo_delta_ppm']:.3f} +/- {row['nnlo_mc_error_ppm']:.3f} ppm")
    print(a.output)


if __name__ == '__main__':
    main()
