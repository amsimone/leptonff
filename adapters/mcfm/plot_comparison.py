#!/usr/bin/env python3
"""Standalone W pT/|eta| and Z mass comparison figures and machine-readable bins."""
import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault('XDG_CACHE_HOME', '/private/tmp/leptonff-cache')
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/leptonff-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

LABELS = {'photonic': ('Analytic photonic NLO + LL', 'C0', '-'),
          'photonic_nll': ('Analytic photonic NLL (unitary)', 'C4', '-.'),
          'leptons': ('Analytic leptonic NLL', 'C2', '--'),
          'qed_hvp': ('Analytic QED + HVP (pointlike)', 'C3', ':')}
TITLES = {'pt': r'$p_{T\ell}$ [GeV]', 'abseta': r'$|\eta_\ell|$',
          'mass': r'$m_{\ell\ell}$ [GeV]'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('report', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    report = json.loads(a.report.read_text())
    a.output.mkdir(parents=True, exist_ok=True)
    models = report['summary']
    primary_model = 'photonic' if 'photonic' in models else 'photonic_nll'
    for process in ['wplus', 'wminus', 'z']:
        observables = ['mass'] if process == 'z' else ['pt', 'abseta']
        for normalised in [False, True]:
            fig, axes = plt.subplots(2, len(observables), figsize=(6*len(observables), 7),
                                     gridspec_kw={'height_ratios': [2.5, 1]},
                                     sharex='col', squeeze=False)
            for col, obs in enumerate(observables):
                primary = models[primary_model][process]['spectra'][obs]['bins']
                finite = [r for r in primary if r['upper'] is not None]
                x = np.array([(r['lower']+r['upper'])/2 for r in finite])
                dx = np.array([(r['upper']-r['lower'])/2 for r in finite])
                widths = 2*dx
                field = 'normalised_fractions' if normalised else 'integrals_pb'
                top, bottom = axes[:, col]
                for component, label, colour in [('born', 'Before FSR', '0.55'),
                                                  ('photos', 'PHOTOS exponentiated', 'C1')]:
                    values = np.array([r[field][component] for r in finite])/widths[:, None]
                    top.errorbar(x, values[:, 0], yerr=values[:, 1], xerr=dx,
                                 fmt='o-', ms=3, lw=1, color=colour, label=label)
                for model, (label, colour, linestyle) in LABELS.items():
                    if model not in models:
                        continue
                    rows = [r for r in models[model][process]['spectra'][obs]['bins'] if r['upper'] is not None]
                    values = np.array([r[field]['match'] for r in rows])/widths[:, None]
                    top.errorbar(x, values[:, 0], yerr=values[:, 1], xerr=dx,
                                 fmt='.', ms=4, lw=1, color=colour, ls=linestyle, label=label)
                if normalised:
                    bottom.axhline(0, color='0.5', lw=.8)
                    bottom.set_ylabel('PHOTOS − analytic\nshape [% / unit]')
                    top.set_ylabel(r'$(1/\sigma_{\rm fid})\,d\sigma/dx$')
                else:
                    bottom.axhline(1, color='0.5', lw=.8)
                    bottom.set_ylabel('PHOTOS / analytic')
                    top.set_ylabel('dσ/dx [pb / unit]')
                for model, (_, colour, _) in LABELS.items():
                    if model not in models:
                        continue
                    rows = [r for r in models[model][process]['spectra'][obs]['bins'] if r['upper'] is not None]
                    marker = {'photonic': 'o', 'photonic_nll': 'D', 'leptons': 's', 'qed_hvp': '^'}[model]
                    if normalised:
                        delta = np.array([r['normalised_photos_minus_match'] for r in rows])*100/widths[:, None]
                        bottom.errorbar(x, delta[:, 0], yerr=delta[:, 1], xerr=dx, fmt=marker, ms=3, color=colour)
                    else:
                        valid = [(k, r['photos_over_match']) for k, r in enumerate(rows) if r['photos_over_match'] is not None]
                        bottom.errorbar([x[k] for k, _ in valid], [r[0] for _, r in valid],
                                        yerr=[r[1] for _, r in valid], xerr=[dx[k] for k, _ in valid],
                                        fmt=marker, ms=3, color=colour)
                bottom.set_xlabel(TITLES[obs])
                for axis in [top, bottom]:
                    axis.grid(alpha=.2)
                if len(finite) != len(primary):
                    overflow = primary[-1]['integrals_pb']['match']
                    top.text(.97, .97, f"Overflow > {primary[-1]['lower']:g} GeV:\n{overflow[0]:.2f} ± {overflow[1]:.2f} pb",
                             transform=top.transAxes, ha='right', va='top', fontsize=8)
            axes[0, 0].legend(fontsize=8, loc='best')
            name = {'wplus': r'$W^+\to\mu^+\nu$', 'wminus': r'$W^-\to\mu^-\bar\nu$',
                    'z': r'$Z/\gamma^*\to\mu^+\mu^-$'}[process]
            samples = models[primary_model][process]['samples_averaged']
            fig.suptitle(name+'; NLO QCD, pp at 13 TeV, bare muons', fontsize=13)
            fig.text(.5, .066, 'pT > 25 GeV, |η| < 2.4'+('; 60 < mℓℓ < 120 GeV' if process == 'z' else ''),
                     ha='center', fontsize=9)
            fig.text(.5, .012, f'{samples} independent seed sample'+('s' if samples != 1 else '')+
                     '; integration errors only.\nPhoton-only LL/NLL: fixed alpha, no pairs; leptonic/HVP curves have different running/pair content.',
                     ha='center', va='bottom', fontsize=8)
            fig.tight_layout(rect=(0, .10, 1, .95))
            stem = process+('_normalised' if normalised else '_comparison')
            for extension in ['png', 'svg', 'pdf']:
                fig.savefig(a.output/f'{stem}.{extension}', dpi=170)
            plt.close(fig)
    with (a.output/'comparison_bins.csv').open('w') as stream:
        writer = csv.writer(stream)
        writer.writerow(['model', 'process', 'observable', 'lower', 'upper', 'component',
                         'integral_pb', 'statistical_error_pb', 'normalised_fraction', 'normalised_error'])
        for model, processes in models.items():
            for process, data in processes.items():
                for obs in (['mass'] if process == 'z' else ['pt', 'abseta']):
                    for row in data['spectra'][obs]['bins']:
                        for component, values in row['integrals_pb'].items():
                            writer.writerow([model, process, obs, row['lower'], row['upper'], component, *values,
                                             *row['normalised_fractions'].get(component, ['', ''])])


if __name__ == '__main__':
    main()
