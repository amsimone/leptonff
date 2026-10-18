#!/usr/bin/env python3
"""Paired PHOTOS ratios, normalised shapes and numerical W/Z controls."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from analyse_spectra import analyse, definitions


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def error(covariance, gradient):
    variance = sum(x*y*covariance[i][j] for i, x in gradient.items()
                   for j, y in gradient.items())
    scale = sum(abs(x*y*covariance[i][j]) for i, x in gradient.items()
                for j, y in gradient.items())
    if variance < -1e-10*max(1e-30, scale):
        raise RuntimeError('Negative propagated variance')
    return math.sqrt(max(0., variance))


def ratio(mean, covariance, numerator, denominator):
    """Delta-method ratio retaining the numerator/denominator covariance."""
    n, d = mean[numerator], mean[denominator]
    if abs(d) <= 3*math.sqrt(max(0., covariance[denominator][denominator])):
        return None  # No useful Gaussian ratio for a denominator compatible with zero.
    return [n/d, error(covariance, {numerator: 1/d, denominator: -n/d**2})]


def aggregate(reports):
    if any(r['labels'] != reports[0]['labels'] for r in reports):
        raise RuntimeError('Cannot average different measurement contracts')
    n = len(reports)
    return ([sum(r['mean_pb'][i] for r in reports)/n for i in range(len(reports[0]['labels']))],
            [[sum(r['covariance_pb2'][i][j] for r in reports)/n**2
              for j in range(len(reports[0]['labels']))]
             for i in range(len(reports[0]['labels']))])


def quadratic_diagnostic(delta, covariance, relative_cut=1e-10):
    """PSD eigenmode chi-square, scaled before resolving exact constraints."""
    import numpy as np
    matrix = np.asarray(covariance, dtype=float)
    vector = np.asarray(delta, dtype=float)
    if not np.isfinite(matrix).all() or not np.isfinite(vector).all():
        raise RuntimeError('Nonfinite diagnostic input')
    scale = float(np.max(np.abs(matrix)))
    if not scale:
        return {'chi2': 0. if not np.any(vector) else None, 'covariance_rank': 0,
                'null_projection_norm': float(np.linalg.norm(vector))}
    matrix = (matrix/scale+matrix.T/scale)/2
    values, modes = np.linalg.eigh(matrix)
    largest = max(0., float(values[-1]))
    if float(values[0]) < -1e-8*max(largest, 1e-30):
        raise RuntimeError('Indefinite diagnostic covariance')
    scaled = vector/math.sqrt(scale)
    # Explicit scalar projections avoid forming an ill-conditioned inverse
    # across exact zero modes (normalisation and forbidden-tail constraints).
    projections = [math.fsum(float(modes[i, j])*float(scaled[i]) for i in range(len(vector)))
                   for j in range(len(vector))]
    retained = [j for j, value in enumerate(values) if value > relative_cut*largest]
    discarded = set(range(len(vector)))-set(retained)
    chi2 = math.fsum(projections[j]**2/float(values[j]) for j in retained)
    return {'chi2': chi2, 'covariance_rank': len(retained),
            'null_projection_norm': math.sqrt(scale*math.fsum(projections[j]**2 for j in discarded))}


def independent_central_runs(paths, baseline, cache_root):
    """Allow a technical QCD subtraction change, retaining the physics setup."""
    import re
    base_meta = json.loads((baseline/'manifest.json').read_text())
    def environment(meta):
        return {k: v for k, v in meta['environment'].items()
                if k not in ['LEPTONFF_HADRONIC_FILE', 'LEPTONFF_FLAVOUR_FILE', 'LEPTONFF_HARD_SAMPLES']}
    def configuration(path):
        section, lines, alpha, seed = '', [], None, None
        for line in (path/'input.ini').read_text().splitlines():
            value = line.strip()
            if value.startswith('['):
                section = value
            if value.startswith(('runstring =', 'rundir =')):
                continue
            if value.startswith('seed ='):
                seed = int(value.split('=', 1)[1]); continue
            if section == '[dipoles]' and value.startswith('alpha ='):
                alpha = float(value.split('=', 1)[1]); continue
            lines.append(line)
        return lines, alpha, seed
    base_config, _, _ = configuration(baseline)
    seeds, alphas, reports, evidence = [], [], [], {}
    for item in paths:
        path = Path(item)
        meta = json.loads((path/'manifest.json').read_text())
        if digest(path/'input.ini') != meta['config_sha256']:
            raise RuntimeError('Changed ensemble configuration: '+str(path))
        if meta.get('hadronic_input') and digest(path/'hadronic.dat') != meta['hadronic_input']['file_sha256']:
            raise RuntimeError('Changed ensemble hadronic input: '+str(path))
        config, alpha, seed = configuration(path)
        if (meta['status'] != 'completed' or environment(meta) != environment(base_meta)
                or meta['requested_inputs'] != base_meta['requested_inputs'] or config != base_config
                or meta.get('hadronic_input') != base_meta.get('hadronic_input')
                or meta['executable_sha256'] != base_meta['executable_sha256']
                or meta['process'] != base_meta['process']):
            raise RuntimeError('Central ensemble changes the physics setup: '+str(path))
        seeds.append(seed); alphas.append(alpha)
        expected = {name: digest(path/name) for name in ['manifest.json', 'iteration_rates.csv', 'event_covariance.csv']}
        cache = cache_root/f'{path.name}_ensemble.json'
        report = json.loads(cache.read_text()) if cache.exists() else None
        if (report is None or report['provenance'] != expected
                or report['analysis_source_sha256'] != digest(Path(__file__).with_name('analyse_spectra.py'))):
            report = analyse(path)
            cache.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
        reports.append(report)
        evidence[path.name] = {'run_directory': str(path.resolve()), 'files': expected,
            'executable_sha256': meta['executable_sha256'], 'source_sha256': meta['source_sha256'],
            'seed': seed, 'dipole_alpha': alpha, 'warnings': meta['integration_warnings']}
    if len(set(seeds)) != len(seeds) or None in seeds or len(set(alphas)) != 1:
        raise RuntimeError('Central ensemble repeats seeds or mixes subtraction parameters')
    return reports, evidence


def rows_for(process, labels, mean, covariance, chunk_covariance=None, seed_covariance=None):
    index = {label: i for i, label in enumerate(labels)}
    total = dict(zip(['born', 'nlo', 'match', 'order2', 'higher', 'psborn', 'photos', 'psdiff'], labels[:8]))
    result = {}
    for observable, edges in definitions(process).items():
        rows = []
        gradients = []
        for k, (lower, upper) in enumerate(zip(edges[:-1], edges[1:])):
            indices = {c: index[f'sp_{observable}{k+1:02d}_{c}']
                       for c in ['born', 'nlo', 'match', 'higher', 'psborn', 'photos', 'psdiff']}
            row = {'lower': lower, 'upper': upper,
                   'integrals_pb': {c: [mean[i], math.sqrt(max(0., covariance[i][i]))]
                                    for c, i in indices.items()},
                   'photos_over_match': ratio(mean, covariance, indices['photos'], indices['match'])}
            if chunk_covariance is not None:
                row['photos_over_match_iteration_diagnostic'] = ratio(mean, chunk_covariance,
                    indices['photos'], indices['match'])
                row['match_iteration_to_event_variance_ratio'] = (
                    chunk_covariance[indices['match']][indices['match']]/covariance[indices['match']][indices['match']]
                    if covariance[indices['match']][indices['match']] else None)
            if seed_covariance is not None:
                row['photos_over_match_seed_diagnostic'] = ratio(mean, seed_covariance,
                    indices['photos'], indices['match'])
            g = {indices['photos']: 1., indices['psborn']: -1.,
                 indices['match']: -1., indices['born']: 1.}
            row['radiation_difference_pb'] = [sum(mean[i]*x for i, x in g.items()), error(covariance, g)]
            shape_gradient = {}
            row['normalised_fractions'] = {}
            for component in ['born', 'nlo', 'match', 'photos']:
                i, t = indices[component], index[total[component]]
                if mean[t] == 0:
                    raise RuntimeError('Zero fiducial denominator')
                g = {i: 1/mean[t], t: -mean[i]/mean[t]**2}
                row['normalised_fractions'][component] = [mean[i]/mean[t], error(covariance, g)]
                if component in ['photos', 'match']:
                    sign = 1 if component == 'photos' else -1
                    for j, x in g.items():
                        shape_gradient[j] = shape_gradient.get(j, 0.)+sign*x
            row['normalised_photos_minus_match'] = [row['normalised_fractions']['photos'][0]
                                                     -row['normalised_fractions']['match'][0],
                                                     error(covariance, shape_gradient)]
            gradients.append(shape_gradient)
            rows.append(row)
        # Include cross-bin correlations in the exported shape difference.
        shape_cov = [[sum(x*y*covariance[i][j] for i, x in g.items() for j, y in h.items())
                      for h in gradients] for g in gradients]
        # Identical tool shapes can cancel to roundoff. Scale the closure
        # tolerance by the component variances before that cancellation.
        shape_scale = max(1e-30, max(row['normalised_fractions'][c][1]**2
                                    for row in rows for c in ['photos', 'match']))
        if max(abs(sum(r)) for r in shape_cov) > 1e-7*shape_scale:
            raise RuntimeError('Shape difference covariance violates the sum constraint')
        result[observable] = {'bins': rows, 'normalised_difference_covariance': shape_cov}
    return result


def chunk_difference(left, right, component='match'):
    """Same-iteration paired diagnostic; NOT an event-product precision bound."""
    def read(path):
        result = {}
        total_labels = ['photos_post'] if component == 'photos' else ['dressed_matched', 'coupled_matched']
        with (path/'iteration_rates.csv').open() as stream:
            for r in csv.reader(stream):
                if ((r[5].startswith('sp_') and r[5].endswith('_'+component))
                        or r[5] in total_labels):
                    result[(tuple(r[:5]), r[5])] = float(r[6])*1e-3
        return result
    l, r = read(left), read(right)
    if l.keys() != r.keys():
        raise RuntimeError('Paired chunks require identical iteration and call counts')
    groups = {}
    for (key, label), value in l.items():
        groups.setdefault((key[:3], label), []).append(r[key, label]-value)
    means, variances = {}, {}
    for (_, label), values in groups.items():
        n = len(values)
        avg = sum(values)/n
        means[label] = means.get(label, 0.)+avg
        variances[label] = variances.get(label, 0.)+sum((x-avg)**2 for x in values)/(n*(n-1))
    return {label: [mean, math.sqrt(variances[label])] for label, mean in means.items()}


def main():
    import numpy as np
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('results', type=Path)
    p.add_argument('--models', nargs='+', default=['photonic', 'leptons', 'qed_hvp'])
    p.add_argument('--cases', nargs='+', default=['nominal', 'seed', 'grid', 'domain', 'cutoff'])
    p.add_argument('--extra-photonic-cases', nargs='*', default=[])
    p.add_argument('--photonic-central', default='nominal',
                   help='Use an included higher-sampling case instead of nominal for the seed average')
    p.add_argument('--photonic-seed', default='seed',
                   help='Independent companion to the photonic central case')
    p.add_argument('--domain-results', type=Path,
                   help='Corrected photonic domain/tail runs from the adapter control-routing fix')
    p.add_argument('--adapter-equivalence', type=Path,
                   help='Required bitwise default-regression evidence when domain controls use a different adapter')
    p.add_argument('--central-runs', type=Path,
                   help='Optional model/process mapping to an independent central seed ensemble')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    checks, flags, summaries, provenance = [], [], {}, {}
    equivalence = json.loads(a.adapter_equivalence.read_text()) if a.adapter_equivalence else None
    overrides = json.loads(a.central_runs.read_text()) if a.central_runs else {}
    if a.domain_results and (not equivalence or not equivalence['passed']):
        raise RuntimeError('Corrected adapter controls require passing default-regression evidence')
    def check(label, passed, **details):
        checks.append({'label': label, 'passed': bool(passed), **details})
    for model in a.models:
        summaries[model] = {}
        cases = list(dict.fromkeys(a.cases+(a.extra_photonic_cases if model == 'photonic' else [])))
        for process in ['wplus', 'wminus', 'z']:
            paths = {case: a.results/f'{model}_{process}_{case}' for case in cases}
            if model == 'photonic' and a.domain_results:
                for case in ['domain', 'tail']:
                    if case in paths:
                        paths[case] = a.domain_results/f'{model}_{process}_{case}'
            reports, metas = {}, {}
            for case, path in paths.items():
                meta = json.loads((path/'manifest.json').read_text())
                if digest(path/'input.ini') != meta['config_sha256']:
                    raise RuntimeError('Changed run configuration: '+str(path))
                if meta.get('hadronic_input') and digest(path/'hadronic.dat') != meta['hadronic_input']['file_sha256']:
                    raise RuntimeError('Changed hadronic input: '+str(path))
                metas[case] = meta
                env = meta['environment']
                if model in ['photonic', 'photonic_nll'] and case in ['domain', 'tail']:
                    if model == 'photonic' and a.domain_results and meta['executable_sha256'] != equivalence['new_executable_sha256']:
                        raise RuntimeError('Unexpected adapter in corrected domain control')
                    requested_grid = {'bins': int(env['LEPTONFF_GRID_BINS']),
                                      'log_span': float(env['LEPTONFF_LOG_SPAN'])}
                    # Older launch scripts did not record resolved_grid. The
                    # corrected bridge always logs the applied values.
                    import re
                    applied = re.findall(r'grid=(\d+), log_span=([^,]+), mu_frag=', (path/'run.log').read_text())
                    if not applied or requested_grid != {'bins': int(applied[-1][0]), 'log_span': float(applied[-1][1])}:
                        raise RuntimeError('Requested photonic domain was not applied')
                expected_photonic = model in ['photonic', 'photonic_nll']
                if (env['LEPTONFF_PHOTONIC_LL'] == '1') != expected_photonic:
                    raise RuntimeError('Wrong theory branch')
                if (env.get('LEPTONFF_PHOTONIC_NLL', '0') == '1') != (model == 'photonic_nll'):
                    raise RuntimeError('Wrong photonic logarithmic accuracy')
                if model == 'photonic_nll' and not meta.get('photonic_nll', {}).get('diagnostics'):
                    raise RuntimeError('Missing applied photonic NLL convention and diagnostics')
                if model == 'leptons' and env['LEPTONFF_FLAVOURS'] != 'leptons':
                    raise RuntimeError('Wrong lepton flavour prescription')
                if model == 'qed_hvp' and (env['LEPTONFF_FLAVOURS'] != 'hadronic-pointlike'
                                           or env['LEPTONFF_QCD_OFF'] != '1'):
                    raise RuntimeError('Expected pure-QED HVP prescription without strong FF evolution')
                cache = a.results/f'{model}_{process}_{case}.json'
                expected = {name: digest(path/name) for name in ['manifest.json', 'iteration_rates.csv', 'event_covariance.csv']}
                report = json.loads(cache.read_text()) if cache.exists() else None
                if (report is None or report['provenance'] != expected
                        or report['analysis_source_sha256'] != digest(Path(__file__).with_name('analyse_spectra.py'))):
                    report = analyse(path)
                    cache.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
                reports[case] = report
                check(f'{model} {process} {case}: grouped measurements and covariance', True,
                      bin_rate_residual_pb=report['maximum_bin_rate_or_overlap_residual_pb'],
                      covariance_closure_pb2=report['maximum_covariance_bin_rate_residual_pb2'])
                if meta['integration_warnings']:
                    flags.append({'model': model, 'process': process, 'case': case,
                                  'warnings': meta['integration_warnings'],
                                  'maximum_chisq_per_iteration': meta['maximum_chisq_per_iteration']})
                provenance[f'{model}_{process}_{case}'] = {'files': expected,
                    'run_directory': str(path.resolve()), 'run_log_sha256': digest(path/'run.log'),
                    'executable_sha256': meta['executable_sha256'], 'source_sha256': meta['source_sha256'],
                    'requested_inputs': meta['requested_inputs'], 'environment': env}
            executables = {m['executable_sha256'] for m in metas.values()}
            if len(executables) != 1:
                if (model != 'photonic' or not equivalence or not equivalence['passed']
                        or executables != {equivalence['old_executable_sha256'], equivalence['new_executable_sha256']}):
                    raise RuntimeError('Mixed executables without verified default equivalence')
            central_case = a.photonic_central if model == 'photonic' else 'nominal'
            companion = a.photonic_seed if model == 'photonic' else 'seed'
            if model == 'photonic' and a.photonic_seed != 'seed' and companion not in reports:
                raise RuntimeError('The requested photonic companion must be included')
            central_cases = list(dict.fromkeys(c for c in [central_case, companion] if c in reports))
            if central_case not in reports:
                raise RuntimeError('The requested central case must be included')
            pooled = [reports[c] for c in central_cases]
            if len(pooled) == 2:
                nominal_env, seed_env = metas[central_case]['environment'], metas[companion]['environment']
                def input_environment(env):
                    return {k: v for k, v in env.items()
                            if k not in ['LEPTONFF_HADRONIC_FILE', 'LEPTONFF_FLAVOUR_FILE', 'LEPTONFF_HARD_SAMPLES']}
                if (input_environment(nominal_env) != input_environment(seed_env)
                        or metas[central_case]['requested_inputs'] != metas[companion]['requested_inputs']
                        or metas[central_case].get('hadronic_input') != metas[companion].get('hadronic_input')):
                    raise RuntimeError('Different physics or budget in independent-seed samples')
                # Validate that the input seed differs; a directory name is not evidence of independence.
                configs = [(paths[c]/'input.ini').read_text() for c in central_cases]
                if configs[0].split('seed = ')[1].splitlines()[0] == configs[1].split('seed = ')[1].splitlines()[0]:
                    raise RuntimeError('Independent-seed average repeats the same seed')
            override_paths = overrides.get(model, {}).get(process)
            ensemble_info = None
            if override_paths:
                pooled, ensemble_info = independent_central_runs(override_paths, paths['nominal'], a.results)
                central_cases = list(ensemble_info)
                for name, info in ensemble_info.items():
                    check(name+': independent ensemble physics, grouped measurements and covariance', True)
                    provenance['ensemble_'+name] = info
                    if info['warnings']:
                        flags.append({'model': model, 'process': process, 'case': name, 'warnings': info['warnings']})
            means, cov = aggregate(pooled)
            seed_cov = ([[sum((r['mean_pb'][i]-means[i])*(r['mean_pb'][j]-means[j]) for r in pooled)
                          /(len(pooled)*(len(pooled)-1)) for j in range(len(means))]
                         for i in range(len(means))] if len(pooled) > 1 else None)
            chunk_cov = [[sum(r['chunk_covariance_diagnostic_pb2'][i][j] for r in pooled)/len(pooled)**2
                          for j in range(len(means))] for i in range(len(means))]
            labels = pooled[0]['labels']
            central = rows_for(process, labels, means, cov, chunk_cov, seed_cov)
            for data in central.values():
                delta = np.array([row['normalised_photos_minus_match'][0] for row in data['bins']])
                matrix = np.array(data['normalised_difference_covariance'])
                data['normalised_shape_comparison'] = {
                    **quadratic_diagnostic(delta, matrix),
                    'scope': 'PHOTOS minus analytic normalised shape; paired event covariance; algorithm configurations differ'}
            totals = {label: [means[i], math.sqrt(max(0., cov[i][i]))] for i, label in enumerate(labels[:8])}
            sensitivities, total_sensitivities = {}, {}
            for case in cases:
                if case in ['nominal', 'seed']:
                    continue
                reference = 'seed' if case == 'angular_seed' else 'nominal'
                component = 'photos' if case == 'photos_cutoff' else 'match'
                rate_index = 6 if component == 'photos' else 2
                paired = chunk_difference(paths[reference], paths[case], component)
                total_sensitivities[case] = {'reference_case': reference,
                    'component': component,
                    'shift_pb': reports[case]['mean_pb'][rate_index]-reports[reference]['mean_pb'][rate_index],
                    'separate_run_quadrature_error_pb': math.sqrt(max(0., reports[case]['covariance_pb2'][rate_index][rate_index]
                                                                         +reports[reference]['covariance_pb2'][rate_index][rate_index])),
                    'paired_iteration_diagnostic_pb': paired[labels[rate_index]]}
                controls = {}
                for obs, edges in definitions(process).items():
                    controls[obs] = []
                    for k in range(len(edges)-1):
                        label = f'sp_{obs}{k+1:02d}_{component}'
                        i = labels.index(label)
                        shift = reports[case]['mean_pb'][i]-reports[reference]['mean_pb'][i]
                        sigma = math.hypot(math.sqrt(max(0., reports[case]['covariance_pb2'][i][i])),
                                           math.sqrt(max(0., reports[reference]['covariance_pb2'][i][i])))
                        controls[obs].append({'lower': edges[k], 'upper': edges[k+1], 'shift_pb': shift,
                            'reference_case': reference,
                            'component': component,
                            'nominal_born_pb': reports[reference]['mean_pb'][labels.index(f'sp_{obs}{k+1:02d}_born')],
                            'reference_component_pb': reports[reference]['mean_pb'][i],
                            'separate_run_quadrature_error_pb': sigma,
                            'paired_iteration_diagnostic_pb': paired[label],
                            'shift_over_nominal_born': shift/reports[reference]['mean_pb'][labels.index(f'sp_{obs}{k+1:02d}_born')]
                            if reports[reference]['mean_pb'][labels.index(f'sp_{obs}{k+1:02d}_born')] else None})
                sensitivities[case] = controls
                # These numerical changes must preserve the Born measurement exactly.
                if case in ['grid', 'domain', 'tail', 'cutoff', 'angular', 'angular_seed', 'photos_cutoff']:
                    residual = max(abs(reports[case]['mean_pb'][i]-reports[reference]['mean_pb'][i])
                                   for i, label in enumerate(labels) if label.endswith('_born') or label == 'ff_born')
                    check(f'{model} {process} {case}: Born stream invariant', residual < 1e-9, residual_pb=residual)
            seed_diagnostic = {}
            if len(pooled) >= 2:
                for obs in central:
                    indices = [labels.index(f'sp_{obs}{k+1:02d}_match') for k in range(len(central[obs]['bins']))]
                    delta = np.array([r['mean_pb'][i]-pooled[0]['mean_pb'][i] for r in pooled[1:] for i in indices])
                    matrix = np.array([[pooled[0]['covariance_pb2'][i][j]
                                        +(pooled[r+1]['covariance_pb2'][i][j] if r == s else 0.)
                                        for s in range(len(pooled)-1) for j in indices]
                                       for r in range(len(pooled)-1) for i in indices])
                    seed_diagnostic[obs] = {'cases': central_cases,
                                            **quadratic_diagnostic(delta, matrix),
                                            'difference_pb': delta.tolist(),
                                            'scope': 'Independent-seed diagnostic; Gaussian covariance may miss rare tails'}
            summaries[model][process] = {'samples_averaged': len(pooled), 'central_cases': central_cases,
                'total_rates_pb': totals,
                'spectra': central, 'numerical_sensitivities': sensitivities,
                'total_sensitivities': total_sensitivities, 'seed_diagnostic': seed_diagnostic,
                'central_ensemble': ensemble_info,
                'independent_seed_covariance_diagnostic_pb2': seed_cov,
                'coupling': metas['nominal'].get('resolved_flavour_coupling'),
                'hadronic_input': metas['nominal'].get('hadronic_input')}
    if 'photonic' in summaries:
        for model in summaries:
            if model == 'photonic':
                continue
            for process in summaries[model]:
                left, right = summaries['photonic'][process], summaries[model][process]
                residual = max(abs(l['integrals_pb'][component][0]-r['integrals_pb'][component][0])
                    for obs in left['spectra']
                    for l, r in zip(left['spectra'][obs]['bins'], right['spectra'][obs]['bins'])
                    for component in ['born', 'psborn', 'photos'])
                check(f'{model} {process}: common production Born and PHOTOS reference', residual < 1e-7,
                      maximum_residual_pb=residual)
    failed = [c for c in checks if not c['passed']]
    report = {'required_checks': len(checks), 'failed_checks': len(failed), 'checks': checks,
              'host_warning_flags': flags, 'summary': summaries, 'provenance': provenance,
              'analysis_sha256': digest(Path(__file__)),
              'central_runs_sha256': digest(a.central_runs) if a.central_runs else None,
              'adapter_equivalence': {'file_sha256': digest(a.adapter_equivalence), 'evidence': equivalence}
              if a.adapter_equivalence else None,
              'limitations': ['Integration errors are statistical, not complete theory uncertainties',
                  'Photonic LL and unitary-conditioned photonic NLL share fixed coupling and no-pair content with PHOTOS',
                  'PHOTOS W ME correction on, Z ME correction off; leading-power analytic vs massive PHOTOS kinematics',
                  'Leptons means e/mu/tau only; qed_hvp adds the explicit pointlike hadronic input prescription with QCD FF evolution off',
                  'Cross-run refinement errors use separate-run quadrature; paired ten-iteration errors are diagnostics only',
                  'NNLO PDF input is shared by the NLO QCD matrix elements; no PDF uncertainty claim',
                  'Fixed-order QCD recoil, inclusive tagged/pair definitions and absence of full NLO EW limit comparison with precision data']}
    a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'required_checks': len(checks), 'failed_checks': len(failed),
                      'host_warning_flags': len(flags)}, indent=2))
    if failed:
        raise SystemExit('Comparison consistency checks failed')


if __name__ == '__main__':
    main()
