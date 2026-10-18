#!/usr/bin/env python3
"""Show photonic NLL on a common angular-inclusive DYTurbo mass spectrum.

This uses a bin-integrated production / midpoint decay migration convolution.
Supply a second production grid to quantify its mass-discretisation dependence.
It is restricted to inclusive neutral-current mass: no laboratory muon cuts.
"""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import sha256, validate_result, MODES


def convolve(production, executable, out, grid):
    r = json.loads(production.read_text())
    validate_result(r, MODES['nlo-nll'])
    if r['process'] != 3 or r['fixed_order'] or r['qcd_order'] != 1 or r['fiducial']:
        raise ValueError('Requires angular-inclusive Z NLO+NLL QCD production')
    bounds = [b['mass'] for b in r['bins']]
    if bounds[0][0] != 66 or bounds[-1][1] != 116:
        raise ValueError('This demonstration requires parent masses 66..116 GeV')
    if any(a[1] != b[0] for a, b in zip(bounds, bounds[1:])):
        raise ValueError('Production mass bins must be contiguous')
    table = out.with_suffix('.production.tsv')
    table.write_text(''.join(f"{sum(b['mass'])/2:.17g} {b['value']:.17g} {b['error']:.17g}\n" for b in r['bins']))
    subprocess.run([str(executable.resolve()), str(table), str(out), str(grid)], check=True)
    data = np.loadtxt(out)
    if not np.all(np.isfinite(data)):
        raise ValueError('Non-finite mass response')
    total = sum(b['value'] for b in r['bins'])
    if not np.allclose(data[:, 2:5].sum(axis=0), total, rtol=1e-11):
        raise ValueError('Inclusive decay migration failed to preserve normalization')
    return data, r


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--production', type=Path, required=True)
    p.add_argument('--coarse-production', type=Path)
    p.add_argument('--coarse-spectrum', type=Path, help='Previously convolved coarse production at the refined FF grid')
    p.add_argument('--decay-executable', type=Path, default=Path('build-dyturbo-repro/leptonff_dyturbo_mass'))
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--grid-bins', type=int, default=8192)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    data, production = convolve(a.production, a.decay_executable, a.output / 'spectrum.tsv', a.grid_bins)
    fine, _ = convolve(a.production, a.decay_executable, a.output / 'ff-refined.tsv', 2*a.grid_bins)
    controls = {}
    visible = data[:, 2] > .001 * data[:, 2].max()
    controls['ff_refinement_max_change_over_born'] = float(np.max(np.abs(fine[visible,4]-data[visible,4])/data[visible,2]))
    peak_mask=(fine[:,0]>=90)&(fine[:,1]<=92)
    controls['peak_ff_refinement_change_over_born'] = float(abs(fine[peak_mask,4].sum()-data[peak_mask,4].sum())/fine[peak_mask,2].sum())
    if a.coarse_production or a.coarse_spectrum:
        coarse = np.loadtxt(a.coarse_spectrum) if a.coarse_spectrum else convolve(a.coarse_production, a.decay_executable, a.output / 'production-coarse.tsv', 2*a.grid_bins)[0]
        if coarse.shape != fine.shape or not np.array_equal(coarse[:,:2],fine[:,:2]):
            raise ValueError('Coarse response uses incompatible observed bins')
        # Compare migration ratios to suppress independent production-quadrature
        # shifts. Report the raw spectra too, since this is not a precision fit.
        controls['production_refinement_max_ratio_change'] = float(np.max(np.abs(fine[visible,4]/fine[visible,2]-coarse[visible,4]/coarse[visible,2])))
        controls['peak_production_refinement_ratio_change'] = float(abs(fine[peak_mask,4].sum()/fine[peak_mask,2].sum()-coarse[peak_mask,4].sum()/coarse[peak_mask,2].sum()))
    fig, (ax, ratio, higher) = plt.subplots(3,1,figsize=(9,9),sharex=True,gridspec_kw={'height_ratios':[3,1.2,1.2]})
    x = np.r_[fine[1:,0],fine[-1,1]]
    colors = ['black','#e07a18','#1764ab']
    labels = ['QED off',r'QED through $O(\alpha)$','Matched photonic NLL']
    for j in range(3):
        ax.stairs(fine[1:,j+2]/1000,x,label=labels[j],color=colors[j],linewidth=1.8)
    ax.set_yscale('log');ax.set_ylim(.2,None);ax.set_ylabel(r'$d\sigma/dm_{\mu\mu}$ [pb/GeV]')
    ax.set_title(r'DYTurbo NLO+NLL QCD + leptonff decay, $pp$ at 13 TeV')
    ax.text(.025,.04,'Inclusive muon angles; parent mass 66–116 GeV\nNo muon pT/η cuts; fixed α photonic model',transform=ax.transAxes,fontsize=10)
    ax.legend(frameon=False)
    central=(fine[:,0]+fine[:,1])/2
    mask=(fine[:,2]>0)&(central>66)
    for j in (1,2):ratio.plot(central[mask],fine[mask,j+2]/fine[mask,2],label=labels[j],color=colors[j])
    ratio.axhline(1,color='gray',lw=.8);ratio.set_ylabel('QED / QED off')
    selected=fine[:,3]>0
    higher.plot(central[selected],fine[selected,4]/fine[selected,3],color=colors[2])
    higher.axhline(1,color='gray',lw=.8);higher.set_ylabel(r'NLL / $O(\alpha)$')
    higher.set_xlabel(r'Bare dimuon mass $m_{\mu\mu}$ [GeV]')
    higher.set_xlim(50,116);ratio.grid(alpha=.2);higher.grid(alpha=.2);fig.tight_layout()
    fig.savefig(a.output/'mass_comparison.png',dpi=160);fig.savefig(a.output/'mass_comparison.svg');plt.close(fig)
    peak=(fine[:,0]>=90)&(fine[:,1]<=92)
    window=(fine[:,0]>=66)&(fine[:,1]<=116)
    def region(mask):
        v=fine[mask,2:5].sum(axis=0)
        return dict(zip(['qed_off_pb','order_alpha_pb','photonic_nll_pb'],map(float,v/1000)))
    report=dict(scope='Inclusive Z mass proof of concept; bin-integrated QCD production with midpoint decay migration, refined in both grids',
        qed='Fixed alpha=1/137.035999084, photonic unitary NLL, finite O(alpha) decay matching, no real pairs or running alpha',
        qcd='DYTurbo NLO+NLL, shared native production for all three curves',
        provenance=dict(production_sha256=sha256(a.production),decay_executable_sha256=sha256(a.decay_executable),
          physics=production['physics'],pdf=production['pdf'],member=production['pdf_member']),
        parent_mass_range_GeV=[66,116],observed_bin_width_GeV=1,
        alpha_reference='fixed',fragmentation_scale_GeV=91.1876,ff_bins=2*a.grid_bins,
        production_bins=len(production['bins']),inclusive_normalization_pb=float(fine[:,2].sum()/1000),
        peak_90_92=region(peak),window_66_116=region(window),controls=controls,
        limitations=['No laboratory muon cuts and no W angular coupling',
          'Production outside 66..116 excluded: no feed-in from larger parent masses',
          'Midpoint production/decay convolution is an approximation; supplied grids quantify refinement',
          'QCD errors and numerical-control differences are not complete theory uncertainties'])
    (a.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
