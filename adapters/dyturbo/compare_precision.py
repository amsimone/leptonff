#!/usr/bin/env python3
"""Fine-bin inclusive mass comparison, reusing PHOTOS transfer moments."""
import argparse
import json
from pathlib import Path
import numpy as np
from run import sha256,validate_result,MODES


def density_map(weights, width):
    """Linear conservative quartic map: five neighbouring bin integrals."""
    n = len(weights)
    result = np.zeros((n, 5, n))
    for j in range(n):
        first = min(n-5, max(0, j-2))
        centers = np.arange(first, first+5)-j
        matrix = np.array([((centers+.5)**(k+1)-(centers-.5)**(k+1))/(k+1)
                           for k in range(5)]).T
        result[j, :, first:first+5] = np.linalg.inv(matrix)/width
    return result


def unpack_covariance(rows):
    out = np.zeros((len(rows), 5, 5))
    k = 5
    for a in range(5):
        for b in range(a, 5):
            out[:, a, b] = out[:, b, a] = rows[:, k]
            k += 1
    return out


def region_covariance(t, v, trials, region, before_in_region):
    selected = ~region if before_in_region else region
    mu = t[region].sum(axis=0)
    second = (v[selected]*(trials-1)
              + np.einsum('bi,bj->bij', t[selected], t[selected])).sum(axis=0)
    return (second-np.outer(mu,mu))/(trials-1)


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--production', required=True, type=Path)
    p.add_argument('--analytic', required=True, type=Path)
    p.add_argument('--photos', required=True, nargs='+', type=Path)
    p.add_argument('--control', nargs='*', default=[], type=Path)
    p.add_argument('--photos-control', type=Path, help='Comparison directory made with looser PHOTOS settings')
    p.add_argument('--campaign', type=Path, help='Generation manifest to include in provenance')
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    production = json.loads(a.production.read_text())
    analytic = json.loads(a.analytic.read_text())
    photos = [json.loads(path.read_text()) for path in a.photos]
    validate_result(production,MODES['nlo-nll'])
    if production['process']!=3 or production['physics']['sqrts']!=13000:
        raise ValueError('This figure requires inclusive neutral-current production at 13 TeV')
    if any(s['schema']!='leptonff-mass-precision-1' for s in [analytic]+photos):
        raise ValueError('Unsupported mass response schema')
    if not all(np.isclose(s['alpha'],analytic['alpha'],rtol=1e-13,atol=0) for s in photos):
        raise ValueError('Analytic and PHOTOS decays use different alpha inputs')
    config = lambda s: (s.get('photos_version','3.56'),s['alpha'],s.get('photos_cutoff',1e-5),
                        s.get('photos_exponentiation_precision',1e-4),s['photos_full_me'])
    if len(set(map(config,photos)))!=1:
        raise ValueError('PHOTOS samples use different physics/numerical configurations')
    if len({s['seed'] for s in photos})!=len(photos):
        raise ValueError('Independent PHOTOS seeds are required')
    mean = np.array([b['mean'] for b in analytic['bins']])
    cov = np.array([b['covariance'] for b in analytic['bins']]).reshape(-1, 3, 3)
    lower = np.array([b['lower'] for b in analytic['bins']])
    upper = np.array([b['upper'] for b in analytic['bins']])
    source = production['bins']
    weights = np.array([b['value'] for b in source])
    errors = np.array([b['error'] for b in source])
    old = photos[0]['transfer_moments']
    aggregate = np.array([[int(b['mass'][0] >= s['lower'] and b['mass'][1] <= s['upper'])
                          for b in source] for s in old])
    if not np.all(aggregate.sum(axis=0) == 1):
        raise ValueError('Production and PHOTOS parent grids do not align')
    if not all(s['upper']-s['lower'] == .25 for s in old):
        raise ValueError('This transfer representation requires 0.25 GeV parent strata')
    mapping = density_map(aggregate@weights, .25)@aggregate
    born_map = np.array([[int(b['mass'][0] >= lo and b['mass'][1] <= hi)
                          for b in source] for lo, hi in zip(lower, upper)])
    if not np.allclose(born_map@weights, mean[:, 0], rtol=1e-10, atol=1e-6):
        raise ValueError('Analytic Born spectrum does not match native production')
    estimates, variances, response_maps = [], [], []
    peak_values, peak_variances = [], []
    peak = (lower >= 90) & (upper <= 92)
    for sample in photos:
        if sample['alpha'] != photos[0]['alpha'] or sample['photos_full_me']:
            raise ValueError('Inconsistent PHOTOS configurations')
        on = mean[:, 0].copy(); variance = np.zeros(len(mean)); response_map = born_map.astype(float).copy()
        peak_value = mean[peak, 0].sum(); peak_variance = 0.
        for j, s in enumerate(sample['transfer_moments']):
            if (s['lower'], s['upper']) != (old[j]['lower'], old[j]['upper']):
                raise ValueError('PHOTOS parent strata differ')
            r = np.array(s['response']); t = r[:, :5]; v = unpack_covariance(r)
            coeff_map = .25*mapping[j]
            c = coeff_map@weights
            on += t@c
            variance += np.einsum('i,bij,j->b', c, v, c)
            response_map += t@coeff_map
            # One parent stratum lies in exactly one observed Born bin.
            # Off-diagonal products with that bin have the opposite sign.
            # Recover a region covariance from diagonal transfer products.
            before_in_peak = 90 <= s['lower'] and s['upper'] <= 92
            mu = t[peak].sum(axis=0)
            region_cov = region_covariance(t,v,s['trials'],peak,before_in_peak)
            peak_value += mu@c
            peak_variance += c@region_cov@c
        estimates.append(on);variances.append(variance);response_maps.append(response_map)
        peak_values.append(peak_value);peak_variances.append(peak_variance)
    estimates, variances = np.array(estimates), np.array(variances)
    on = estimates.mean(axis=0); mc_var = variances.sum(axis=0)/len(photos)**2
    photos_map = np.mean(response_maps, axis=0)
    valid = mean[:, 0] > 0
    nll_ratio = np.divide(mean[:, 2], mean[:, 0], out=np.zeros(len(mean)), where=valid)
    photos_ratio = np.divide(on, mean[:, 0], out=np.zeros(len(mean)), where=valid)
    input_derivative = photos_map-photos_ratio[:, None]*born_map
    photos_input_var = np.sum((input_derivative*errors)**2, axis=1)
    nll_var = cov[:, 2, 2]+nll_ratio**2*cov[:, 0, 0]-2*nll_ratio*cov[:, 2, 0]
    controls = {}
    for path in a.control:
        control = json.loads(path.read_text())
        other = np.array([b['mean'] for b in control['bins']])
        controls[str(path)] = float(np.max(np.abs(other[valid, 2]/other[valid, 0]-nll_ratio[valid])))
    a.output.mkdir(parents=True, exist_ok=False)
    x = np.r_[lower[1:], upper[-1]]
    center = (lower+upper)/2
    fig, (ax, ratio, difference) = plt.subplots(3, 1, figsize=(10, 10), sharex=True,
        gridspec_kw={'height_ratios': [3, 1.4, 1.4]})
    for values, label, color in ((mean[:, 0], 'QED off', 'black'),
        (mean[:, 1], r'LeptonFF through $O(\alpha)$', '#dc7b19'),
        (mean[:, 2], 'LeptonFF photonic NLL', '#1764ab'), (on, 'PHOTOS exponentiated', '#b52a45')):
        ax.stairs(values[1:]/250, x, label=label, color=color, linewidth=1.3)
    ax.set_yscale('log');ax.set_ylabel(r'$d\sigma/dm_{\mu\mu}$ [pb/GeV]')
    ax.set_title(r'DYTurbo NLO+NLL QCD, bare muons, $pp$ at 13 TeV')
    ax.text(.02, .03, 'Inclusive muon angles; parent mass 66–116 GeV\nFixed α; no pairs; PHOTOS production ME weights off',
            transform=ax.transAxes, fontsize=9)
    ax.legend(frameon=False, fontsize=10)
    ratio.plot(center[valid], nll_ratio[valid], label='LeptonFF', color='#1764ab')
    ratio.errorbar(center[valid], photos_ratio[valid], yerr=np.sqrt(mc_var[valid])/mean[valid, 0],
                   fmt='.', ms=2, color='#b52a45', label='PHOTOS, MC error')
    ratio.axhline(1, color='gray', lw=.7);ratio.set_ylabel('QED on / off');ratio.legend(frameon=False, fontsize=9)
    difference.errorbar(center[valid], 100*(photos_ratio[valid]-nll_ratio[valid]),
                        yerr=100*np.sqrt(mc_var[valid])/mean[valid, 0],fmt='.',ms=2,color='#b52a45')
    difference.axhline(0, color='gray', lw=.7);difference.set_ylabel('PHOTOS − LeptonFF\n[percentage points]')
    difference.set_xlabel(r'Bare dimuon mass $m_{\mu\mu}$ [GeV]');difference.set_xlim(66,116)
    ratio.grid(alpha=.2);difference.grid(alpha=.2);fig.tight_layout()
    for suffix in ('png', 'pdf', 'svg'):fig.savefig(a.output/f'mass_comparison.{suffix}',dpi=180)
    plt.close(fig)
    peak_off=mean[peak, 0].sum();peak_photo=np.mean(peak_values);peak_err=np.sqrt(sum(peak_variances))/len(photos)
    report = {'scope': 'Inclusive mass; common DYTurbo NLO+NLL QCD, normalized photonic decays',
        'observed_bin_width_GeV': .25,'parent_mass_GeV':[66,116],
        'photos_trials': sum(s['photos_trials'] for s in photos),'independent_seeds':[s['seed'] for s in photos],
        'photos_settings':dict(zip(('version','alpha','cutoff','exponentiation_precision','full_me'),config(photos[0]))),
        'peak_90_92_pb': {'off': peak_off/1000,'order_alpha':mean[peak,1].sum()/1000,
                          'photonic_nll':mean[peak,2].sum()/1000,'photos':peak_photo/1000,'photos_mc_error':peak_err/1000},
        'max_ratio_mc_error':float(np.max(np.sqrt(mc_var[valid])/mean[valid,0])),
        'peak_max_ratio_mc_error':float(np.max(np.sqrt(mc_var[peak])/mean[peak,0])),
        'max_ratio_photos_input_error':float(np.max(np.sqrt(photos_input_var[valid])/mean[valid,0])),
        'max_ratio_nll_input_error':float(np.max(np.sqrt(np.maximum(0,nll_var[valid]))/mean[valid,0])),
        'production_grid_controls':controls,
        'normalization_residual_fb':float(on.sum()-mean[:,0].sum()),
        'seed_pull_rms':float(np.sqrt(np.mean((estimates[:,valid]-on[valid])**2/variances[:,valid]))),
        'provenance':{str(path):sha256(path) for path in [a.production,a.analytic]+a.photos+a.control+[Path(__file__)]},
        'limitations':['No lab muon cuts; no parent-mass feed-in above 116 GeV.',
          'PHOTOS 3.56 universal decay mode: interference/exponentiation on, full production-dependent ME off.',
          'Error bars in figure are PHOTOS MC only. Native quadrature error estimates and refinement shifts are reported separately.',
          'No QCD/QED scale/PDF uncertainties; this is a numerical tool comparison.',
          'PHOTOS transfer density uses conservative 0.25 GeV reconstruction; its source refinement is supplied as a control.']}
    report['peak_corrections_percent']={k:100*(report['peak_90_92_pb'][k]/(peak_off/1000)-1)
                                       for k in ('order_alpha','photonic_nll','photos')}
    report['peak_photos_minus_nll_percentage_points']=100*(peak_photo-mean[peak,2].sum())/peak_off
    report['precision_met_for_photos_mc']=report['max_ratio_mc_error']<.001
    if a.campaign:
        report['provenance'][str(a.campaign)]=sha256(a.campaign)
    if a.photos_control:
        path=a.photos_control/'report.json';other=json.loads(path.read_text())
        spectrum=a.photos_control/'spectrum.tsv';data=np.loadtxt(spectrum)
        if not np.array_equal(data[:,:2],np.c_[lower,upper]) or not np.allclose(data[:,2],mean[:,0],rtol=1e-12):
            raise ValueError('PHOTOS numerical control uses different production/bins')
        pulls=(on[valid]-data[valid,5])/np.sqrt(mc_var[valid]+data[valid,6]**2)
        report['photos_numerical_control']={
            'loose_settings':other.get('photos_settings',dict(version='3.56',cutoff=1e-5,exponentiation_precision=1e-4)),
            'loose_trials':other['photos_trials'],
            'tight_loose_peak_shift_pb':peak_photo/1000-other['peak_90_92_pb']['photos'],
            'tight_loose_peak_error_pb':float(np.hypot(peak_err/1000,other['peak_90_92_pb']['photos_mc_error'])),
            'per_bin_difference_pull_rms':float(np.sqrt(np.mean(pulls**2))),
            'max_abs_per_bin_pull':float(np.max(np.abs(pulls)))}
        report['provenance'].update({str(path):sha256(path),str(spectrum):sha256(spectrum)})
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    np.savetxt(a.output/'spectrum.tsv',np.c_[lower,upper,mean,on,np.sqrt(mc_var),np.sqrt(photos_input_var),np.sqrt(np.maximum(0,nll_var))],
               header='lower upper off order_alpha photonic_nll photos photos_mc_error photos_ratio_input_numerator_error nll_ratio_input_numerator_error (fb)')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':main()
