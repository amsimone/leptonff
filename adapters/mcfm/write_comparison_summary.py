#!/usr/bin/env python3
"""Write the human-readable comparison alongside its plots and covariance."""
import argparse
import json
from pathlib import Path


def pm(pair, digits=2):
    return f'{pair[0]:.{digits}f} ± {pair[1]:.{digits}f}'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('report', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--input-report', type=Path)
    a = p.parse_args()
    report = json.loads(a.report.read_text())
    models = report['summary']
    primary_model = 'photonic' if 'photonic' in models else 'photonic_nll'
    names = {'wplus': 'W+', 'wminus': 'W−', 'z': 'Z/γ*'}
    lines = ['# W/Z bare-muon fragmentation comparison', '',
             'NLO QCD production at pp √s = 13 TeV. Observed bare muons require '
             'pT > 25 GeV and |η| < 2.4; Z pairs also require 60 < mℓℓ < 120 GeV. '
             'There is no W missing-pT or mT cut. W η is folded as |η|.', '',
             'The direct tool comparison is **analytic fixed-alpha photonic NLO + LL '
             'versus exponentiated PHOTOS 3.56**, with the same production inputs, '
             'fixed alpha and no pair production. The leptonic NLL and pointlike '
             'QED+HVP curves show additional content separately. PHOTOS has the W '
             'matrix-element correction enabled and the Z correction disabled.', '',
             '**Errors below are integration statistics. They are not a complete '
             'theory uncertainty or a precision-data certification.**', '',
             '## Fiducial cross sections', '',
             '| Process | Analytic photonic [pb] | PHOTOS [pb] | Paired PHOTOS − analytic [pb] |',
             '| --- | ---: | ---: | ---: |']
    if primary_model == 'photonic_nll':
        lines = [line.replace('photonic NLO + LL', 'unitary-conditioned photonic NLL') for line in lines]
    if 'photonic_nll' in models:
        lines[2:2] = ['The photonic NLL curve uses fixed alpha and no pairs, with a single time-like NLO-kernel insertion and explicit unitary conditioning. Its source-inclusive factor is recorded in each run manifest.', '']
    for process, data in models[primary_model].items():
        values = list(data['total_rates_pb'].values())
        lines.append(f"| {names[process]} | {pm(values[2])} | {pm(values[6])} | {pm(values[7])} |")
    lines += ['', 'The central estimates average independent seeds. Their full event '
              'covariances are added with the square of the averaging weight. The '
              'photonic central cases are '+', '.join(models[primary_model]['wplus']['central_cases'])+
              '; the higher angular-sampling case replaces the original seed-1045 '
              'estimate when selected, and angular_seed similarly replaces the '
              'seed-2089 estimate. Numerical sensitivities compare with the '
              'original nominal case, except angular_seed, which compares with '
              'its own original seed.', '',
              '## Requested differential comparisons', '',
              'Absolute spectra have a paired PHOTOS/analytic ratio. Normalised '
              'spectra have a paired shape difference with cross-bin covariance. '
              'Ratios with a denominator compatible with zero are withheld, while '
              'their signed integrals remain in the data.', '']
    if models[primary_model]['wminus'].get('central_ensemble'):
        n = models[primary_model]['wminus']['samples_averaged']
        lines += [f'The W− central result uses **{n} independent seeds** with '
                  'QCD dipole alpha = 0.01 and 16 hard-radiation samples in each '
                  'fragmentation mode. Original alpha = 1 samples remain '
                  'numerical controls. The smaller subtraction parameter is '
                  'a technical integration choice; the hard process, cuts '
                  'and physical inputs are unchanged and checked. All ensemble '
                  'seeds are retained.', '']
    for process in ['wplus', 'wminus', 'z']:
        lines += [f'### {names[process]}', '',
                  f'![{names[process]} absolute comparison](plots/{process}_comparison.png)', '',
                  f'[Normalised shapes](plots/{process}_normalised.png) · '
                  f'[Vector figure](plots/{process}_comparison.svg) · '
                  f'[PDF figure](plots/{process}_comparison.pdf)', '',
                  '| Observable | Bin | PHOTOS / photonic − 1 [%] |',
                  '| --- | --- | ---: |']
        for observable in (['mass'] if process == 'z' else ['pt', 'abseta']):
            for row in models[primary_model][process]['spectra'][observable]['bins']:
                edge = f"{row['lower']:g}–{row['upper']:g}" if row['upper'] is not None else f">{row['lower']:g}"
                ratio = row['photos_over_match']
                value = pm([100*(ratio[0]-1), 100*ratio[1]]) if ratio else 'withheld (weak denominator)'
                lines.append(f'| {observable} | {edge} | {value} |')
        lines.append('')
    lines += ['## Paired normalised-shape comparison', '',
              'This uses the complete cross-bin covariance of PHOTOS minus the '
              'analytic normalised shape. It quantifies agreement of these '
              'algorithm configurations; it is not an independent physical '
              'accuracy test. The covariance rank accounts for the sum-to-one '
              'constraint. Leptonic/HVP comparisons have different running/pair content.', '',
              '| Model | Process | Observable | χ² / covariance rank |',
              '| --- | --- | --- | ---: |']
    for model, processes in models.items():
        for process, data in processes.items():
            for obs in (['mass'] if process == 'z' else ['pt', 'abseta']):
                value = data['spectra'][obs]['normalised_shape_comparison']
                lines.append(f"| {model} | {names[process]} | {obs} | {value['chi2']:.2f} / {value['covariance_rank']} |")
    lines += ['', '## Numerical controls', '',
              'The following are the largest absolute shifts across the requested '
              'finite bins, relative to the corresponding nominal Born bin. '
              'The accompanying uncertainty is the **paired iteration diagnostic**, '
              'based on ten iterations per QCD contribution, and is not a certified '
              'numerical bound. The JSON also retains separate-run quadrature '
              'errors and every individual bin. The domain control changes both '
              'domain and resolution; the tail control preserves resolution. '
              'The photos_cutoff row measures PHOTOS, while the other rows '
              'measure the analytic matched result.', '',
              '| Model | Process | Control | Largest shift / Born [%] | Paired iteration error [%] |',
              '| --- | --- | --- | ---: | ---: |']
    for model, processes in models.items():
        for process, data in processes.items():
            requested = ['mass'] if process == 'z' else ['pt', 'abseta']
            for case, observables in data['numerical_sensitivities'].items():
                candidates = [r for obs in requested for r in observables[obs]
                              if r['upper'] is not None and r['nominal_born_pb'] != 0]
                worst = max(candidates, key=lambda r: abs(r['shift_over_nominal_born']))
                sigma = 100*worst['paired_iteration_diagnostic_pb'][1]/abs(worst['nominal_born_pb'])
                lines.append(f"| {model} | {names[process]} | {case} | {100*worst['shift_over_nominal_born']:.4f} | {sigma:.4f} |")
    lines += ['', '## Independent-seed diagnostics', '',
              'These compare the selected independent central seed samples, using '
              'their summed event covariance. A large χ² is a reason to investigate '
              'sampling or rare tails; passing algebraic checks does not override it.', '',
              '| Model | Process | Observable | χ² / covariance rank |',
              '| --- | --- | --- | ---: |']
    for model, processes in models.items():
        for process, data in processes.items():
            for obs in (['mass'] if process == 'z' else ['pt', 'abseta']):
                if obs in data['seed_diagnostic']:
                    value = data['seed_diagnostic'][obs]
                    lines.append(f"| {model} | {names[process]} | {obs} | {value['chi2']:.2f} / {value['covariance_rank']} |")
    if a.input_report:
        inputs = json.loads(a.input_report.read_text())
        lines += ['', '## Hadronic-input prescription sensitivity', '',
                  'This auxiliary pure-QED study varies the pointlike matching '
                  'scale from 2 GeV to 1.5/3 GeV and the sourced HVP input by '
                  '±1 input standard deviation. It uses the reduced 0.05 call '
                  'budget and retains '+str(len(inputs['host_warning_flags']))+
                  ' W+ host convergence warnings. The errors here are paired '
                  'iteration diagnostics, not a bound on all hadronic modelling.', '',
                  '| Process | Variation | Total-rate shift [pb] | Paired iteration error [pb] |',
                  '| --- | --- | ---: | ---: |']
        for process, data in inputs['summary']['qed_hvp'].items():
            for case, value in data['total_sensitivities'].items():
                lines.append(f"| {names[process]} | {case} | {value['shift_pb']:.6f} | {value['paired_iteration_diagnostic_pb'][1]:.6f} |")
    lines += ['', '## NLL content and use limits', '',
              '| Model | W+ [pb] | W− [pb] | Z/γ* [pb] |',
              '| --- | ---: | ---: | ---: |']
    for model, processes in models.items():
        values = [pm(list(processes[process]['total_rates_pb'].values())[2]) for process in ['wplus', 'wminus', 'z']]
        lines.append('| '+model+' | '+' | '.join(values)+' |')
    lines += ['', 'The leptonic/HVP models include running/pair content that PHOTOS does not '
              'share. `leptons` contains e/mu/tau; `qed_hvp` uses the explicit '
              'pointlike quark/HVP prescription with strong FF evolution disabled. '
              'Differences between separately run theory models do not have '
              'within-event cross-model covariance and are not assigned a '
              'significance here. The exported radiation difference subtracts '
              'each tool’s own Born map before comparing radiation.', '',
              f"Required consistency checks: **{report['required_checks']-report['failed_checks']}/{report['required_checks']} passed**. "
              f"Host integration warning flags: **{len(report['host_warning_flags'])}**.", '',
              *(['The photonic domain-routing fix is checked against the previous '
                 'adapter: all default rates, spectra and raw/event covariance '
                 'products agree bitwise. Corrected domain controls are explicitly '
                 'identified in provenance; no numerical conclusion is drawn from '
                 'the older photonic runs that ignored the domain setting.', '']
                if report.get('adapter_equivalence') else []),
              'Before comparing with precision data, select the measurement and '
              'match its energy, cuts, prompt/extra-pair definition, line shape '
              'and branching conventions; evaluate PDF/alpha_s and QCD-scale '
              'uncertainties and address recoil resummation near the W Jacobian '
              'edge. Mixed strong/QED fragmentation kernels and fitted hadronic '
              'fragmentation are not prerequisites to the present pure-QED benchmark.', '',
              '[Full machine-readable report](report.json) · '
              '[Bin table](plots/comparison_bins.csv) · '
              '[Accuracy audit](../../docs/qed_accuracy_audit.md)', '']
    a.output.write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
