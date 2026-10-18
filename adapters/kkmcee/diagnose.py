#!/usr/bin/env python3
"""Diagnose fixed-Q residuals with independent coupling and mass controls.

Fits are per observable across independent runs, never across correlated
moments. They diagnose a perturbative expansion, not a theory uncertainty.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from scipy.stats import chi2
from scipy.integrate import quad

from compare import MODES, difference, input_paths, normalized, prediction, records, seeds
from check_mass_radiator import continuum


def fit_powers(factors, values, errors, powers):
    """Fit sum c_p*(alpha/alpha_physical)^p with known independent MC errors.

    No residual-based rescaling of covariance: discrepancies can represent
    omitted powers rather than an incorrectly estimated integration variance.
    """
    f,y,e=map(lambda a:np.asarray(a,dtype=float),(factors,values,errors))
    if not (f.ndim==y.ndim==e.ndim==1 and len(f)==len(y)==len(e) and
            len(f)>=len(powers) and np.all(np.isfinite(f)) and
            np.all(np.isfinite(y)) and np.all(np.isfinite(e)) and np.all(e>0)):
        raise ValueError('Invalid independent fit inputs')
    design=np.column_stack([f**p for p in powers])
    matrix=design/e[:,None]
    u,s,vh=np.linalg.svd(matrix,full_matrices=False)
    if s[-1]<=s[0]*1e-12:raise ValueError('Unresolved coupling fit')
    coefficients=vh.T@((u.T@(y/e))/s)
    covariance=(vh.T/s**2)@vh
    pulls=(y-design@coefficients)/e
    statistic=float(pulls@pulls);dof=len(f)-len(powers)
    return dict(powers=list(powers),coefficients=coefficients.tolist(),
        covariance=covariance.tolist(),errors=np.sqrt(np.diag(covariance)).tolist(),
        chi2=statistic,dof=dof,p_value=float(chi2.sf(statistic,dof)) if dof else None,
        residual_pulls=pulls.tolist())


def support_estimate(d,min_mass):
    """Leading-power O(alpha) estimate only, not a certified tail bound.

    CEEX weights vanish below the upstream MinMassCEEX. Removing that region
    raises normalized N>1 moments; the sign is opposite to the observed deficit.
    The leading-power spectrum is least accurate near the physical threshold.
    """
    lo=(2*d['mass']/d['Q'])**2;hi=(min_mass/d['Q'])**2
    if hi>=1:raise ValueError('CEEX support floor excludes the benchmark')
    L=2*math.log(d['Q']/d['mass'])
    def density(z):return d['alpha']/math.pi*(1+z*z)/(1-z)*(L+math.log(z)-1)
    fraction=quad(density,lo,hi,epsabs=1e-13)[0] if hi>lo else 0.
    return dict(min_mass=min_mass,first_order_missing_fraction_estimate=fraction,
        warning='Leading-power first-order estimate; threshold mass terms and higher orders are not bounded.')


def validate_controls(central, coupling, mass):
    """Only the intended parameter can change; all streams must be disjoint."""
    common=('Q','cutoff','fsr','isr','ifi','weak','pairs','running_alpha','modes',
            'remove_soft_photons','normalization')
    seen=set(seeds(central))
    for d,kind in [(d,'coupling') for d in coupling]+[(mass,'mass')]:
        if any(d.get(k)!=central.get(k) for k in common):
            raise ValueError('Control changes benchmark settings')
        if (d.get('ceex_min_mass',1.)!=central.get('ceex_min_mass',1.) or
                d.get('ceex_active',d['fsr'])!=central.get('ceex_active',central['fsr'])):
            raise ValueError('Control changes CEEX settings or support')
        fixed='mass' if kind=='coupling' else 'alpha'
        if d[fixed]!=central[fixed]:raise ValueError('Control changes an unintended parameter')
        if d.get('provenance',{}).get('defaults_sha256')!=central.get('provenance',{}).get('defaults_sha256'):
            raise ValueError('Control changes KKMC defaults')
        # Old small-alpha controls have a metadata/no-radiation driver update.
        # Their upstream physics build is unchanged; record every binary hash.
        if seen&seeds(d):raise ValueError('Control streams overlap')
        seen.update(seeds(d))
        if [(b['observable'],b['value']) for b in d['observables']]!=[
                (b['observable'],b['value']) for b in central['observables']]:
            raise ValueError('Control observables differ')
    if not central['fsr'] or mass['mass']>=central['mass']:
        raise ValueError('Need an FSR-on small-mass control')
    factors=[d['alpha']/central['alpha'] for d in coupling]
    if len(set(factors+[1.]))!=len(factors)+1:
        raise ValueError('Coupling controls must be distinct from the central point')


def analyze(central,coupling,mass,mass_coupling=()):
    validate_controls(central,coupling,mass)
    scan=sorted([central,*coupling],key=lambda d:d['alpha'])
    seen=set().union(*(seeds(d) for d in [central,*coupling,mass]))
    for d in mass_coupling:
        reference=next(c for c in scan if c['alpha']==d['alpha'])
        validate_controls(reference,[],d)
        if d['mass']!=mass['mass'] or seen&seeds(d):
            raise ValueError('Invalid or overlapping mass-coupling control')
        seen.update(seeds(d))
    rows=[]
    for index,b in enumerate(central['observables']):
        if b['observable']!='moment':continue
        N=b['value'];points=[]
        for d in scan:
            ref=prediction(dict(Q=d['Q'],mass=d['mass'],scale=d['Q'],alpha=d['alpha'],N=N))
            obs=d['observables'][index];aa=d['alpha']/(2*math.pi)
            models={name:normalized(obs,k) for k,name in enumerate(MODES)}
            points.append(dict(alpha_factor=d['alpha']/central['alpha'],trials=d['trials'],
                modes=models,ff=ref,residual=models['ceex2']['value']-ref['nll'],
                mc_error=models['ceex2']['mc_error'],
                ceex2_minus_ceex1=difference(obs,4,3),ceex2_minus_eex2=difference(obs,4,1),
                paired_order_coefficient={k:v/aa**2 for k,v in difference(obs,4,3).items()}))
        factors=[p['alpha_factor'] for p in points]
        fits={}
        for label,powers in [('quadratic',[2]),('quadratic_cubic',[2,3]),
                              ('linear_quadratic_cubic',[1,2,3]),
                              ('linear_through_quartic',[1,2,3,4])]:
            fits[label]=fit_powers(factors,[p['residual'] for p in points],
                                  [p['mc_error'] for p in points],powers)
        # A separate paired-weight test can resolve lower-order differences
        # between the EEX/CEEX schemes much more precisely than absolute ratios.
        fits['ceex2_minus_eex2']=fit_powers(factors,
            [p['ceex2_minus_eex2']['value'] for p in points],
            [p['ceex2_minus_eex2']['mc_error'] for p in points],[1,2,3])
        mass_points=[]
        for d in (central,mass):
            a=d['alpha']/(2*math.pi);L=2*math.log(d['Q']/d['mass'])
            ref=prediction(dict(Q=d['Q'],mass=d['mass'],scale=d['Q'],alpha=d['alpha'],N=N))
            mc=normalized(d['observables'][index],4)
            # Strip the common LL damping before testing a finite alpha^2
            # remainder. At finite coupling, alpha^3*L and higher also contribute.
            factor=a*a*math.exp(2*a*L*continuum(N)['P'])
            mass_points.append(dict(mass=d['mass'],L=L,trials=d['trials'],
                residual=mc['value']-ref['nll'],mc_error=mc['mc_error'],
                reduced_residual=(mc['value']-ref['nll'])/factor,
                reduced_mc_error=mc['mc_error']/factor))
        x,y=mass_points;delta=y['reduced_residual']-x['reduced_residual']
        error=math.hypot(x['reduced_mc_error'],y['reduced_mc_error'])
        delta_L=y['L']-x['L']
        mass_scan=[]
        for d in sorted([mass,*mass_coupling],key=lambda d:d['alpha']):
            physical=next(c for c in scan if c['alpha']==d['alpha'])
            x=normalized(d['observables'][index],4)
            y=normalized(physical['observables'][index],4)
            def residual(record,value):
                return value-prediction(dict(Q=record['Q'],mass=record['mass'],
                    scale=record['Q'],alpha=record['alpha'],N=N))['nll']
            mass_scan.append(dict(alpha_factor=d['alpha']/central['alpha'],
                value=residual(d,x['value'])-residual(physical,y['value']),
                mc_error=math.hypot(x['mc_error'],y['mc_error'])))
        mass_fits={}
        if len(mass_scan)>=3:
            for label,powers in [('quadratic',[2]),('cubic',[3]),('quadratic_cubic',[2,3]),
                                  ('cubic_quartic',[3,4]),('quadratic_through_quartic',[2,3,4])]:
                mass_fits[label]=fit_powers([x['alpha_factor'] for x in mass_scan],
                    [x['value'] for x in mass_scan],[x['mc_error'] for x in mass_scan],powers)
        rows.append(dict(N=N,coupling_points=points,fits=fits,mass_points=mass_points,
            mass_coupling_points=mass_scan,mass_coupling_fits=mass_fits,
            reduced_mass_change=dict(value=delta,mc_error=error,pull=delta/error,
                slope=delta/delta_L,slope_mc_error=error/delta_L)))
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--central',type=Path,nargs='+',required=True)
    p.add_argument('--coupling',type=Path,nargs='+',action='append',required=True,
                   help='Repeat once for each independent coupling point')
    p.add_argument('--small-mass',type=Path,nargs='+',required=True)
    p.add_argument('--mass-coupling',type=Path,nargs='+',action='append',default=[],
                   help='Small-mass controls at the coupling-scan values')
    p.add_argument('--Q',type=float,default=91.1876)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    central=next(d for d in records(a.central) if d['Q']==a.Q)
    coupling=[]
    for folders in a.coupling:
        data=records(folders)
        if len(data)!=1 or data[0]['Q']!=a.Q:raise ValueError('Need one Q per control')
        coupling.append(data[0])
    mass=records(a.small_mass)
    if len(mass)!=1 or mass[0]['Q']!=a.Q:raise ValueError('Need one mass control at the same Q')
    mass_coupling=[]
    for folders in a.mass_coupling:
        data=records(folders)
        if len(data)!=1 or data[0]['Q']!=a.Q:raise ValueError('Need one Q per mass/coupling control')
        mass_coupling.append(data[0])
    rows=analyze(central,coupling,mass[0],mass_coupling)
    folders=[*a.central,*[f for group in a.coupling for f in group],
             *a.small_mass,*[f for group in a.mass_coupling for f in group]]
    backends={}
    for folder in folders:
        path=folder/'backend-manifest.json';backend=json.loads(path.read_text())
        backends[str(path)]=dict(sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            commit=backend['commit'],source_tree_sha256=backend['source_tree_sha256'],
            driver_sha256=backend['driver_sha256'],fortran_flags=backend['fortran_flags'])
    if len({(d['commit'],d['source_tree_sha256'],d['fortran_flags']) for d in backends.values()})!=1:
        raise ValueError('Diagnostic campaigns use different upstream physics builds')
    defaults=(a.central[0]/'defaults.input').read_text()
    min_mass=float(re.search(r'^\s*638\s+(\S+)',defaults,re.M)[1])
    report=dict(schema='leptonff-kkmcee-residual-diagnostics-v1',Q=a.Q,rows=rows,
        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        hypotheses={
            'linear':'a residual proportional to alpha indicates a first-order mismatch',
            'quadratic':'a residual proportional to alpha^2 is consistent with a second-order remainder',
            'mass':'a finite alpha^2 remainder is independent of L after stripping LL damping at leading order',
            'limit':'finite-coupling mass dependence also receives alpha^3*L and higher remainders; one mass lever arm cannot identify a unique kernel'},
        caveats=['Only MC errors enter the fits; polynomial truncation is a separate limitation.',
                 'Fits are per moment across independent coupling points, not a global test across correlated moments.',
                 'CEEX/EEX residual orders are not fragmentation LL/NLL accuracy labels.',
                 'No kernel or matching coefficient is fitted into the production prediction.',
                 'The legacy CEEX denominator excludes dimuon masses below MinMassCEEX; see support estimates.'],
        backend_manifests=backends,
        samples=[dict(Q=d['Q'],mass=d['mass'],alpha=d['alpha'],trials=d['trials'],
                      seeds=sorted(seeds(d)),provenance=d.get('provenance'),
                      ceex_support=support_estimate(d,min_mass))
                 for d in [central,*coupling,mass[0],*mass_coupling]],
        input_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest()
                      for folders in [a.central,*a.coupling,a.small_mass,*a.mass_coupling] for path in input_paths(folders)})
    a.output.mkdir(parents=True,exist_ok=False)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for row,color in zip(rows,('#2563eb','#dc2626','#059669','#9333ea')):
        points=[x for x in row['coupling_points'] if x['alpha_factor']>=.5]
        f=np.array([x['alpha_factor'] for x in points])
        axes[0].errorbar(f,1e6*np.array([x['residual'] for x in points])/f**2,
            yerr=1e6*np.array([x['mc_error'] for x in points])/f**2,
            fmt='o-',color=color,label=f'N = {row["N"]:g}')
        m=row['mass_points']
        axes[1].errorbar([x['L'] for x in m],[x['reduced_residual'] for x in m],
            yerr=[x['reduced_mc_error'] for x in m],fmt='o-',color=color)
    axes[0].set_xscale('log',base=2)
    axes[0].set_xlabel(r'$f=\alpha/\alpha_{\rm physical}$')
    axes[0].set_ylabel(r'$(R_{\rm CEEX2}-R_{\rm NLL})/f^2$ [$10^{-6}$]')
    axes[0].set_title('Coupling dependence; constant = pure quadratic residual')
    axes[0].legend()
    axes[1].set_xlabel(r'$L=\log(Q^2/m^2)$')
    axes[1].set_ylabel(r'$(R_{\rm CEEX2}-R_{\rm NLL})/(a^2 U_N)$')
    axes[1].set_title('Mass dependence; finite-coupling diagnostic')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Fixed-Q photonic residuals; independent runs; MC errors only')
    for suffix in ('png','pdf','svg'):fig.savefig(a.output/('residual-scaling.'+suffix),dpi=180)
    plt.close(fig)
    if mass_coupling:
        fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
        for row,color in zip(rows,('#2563eb','#dc2626','#059669','#9333ea')):
            points=row['mass_coupling_points'];f=np.array([x['alpha_factor'] for x in points])
            ax.errorbar(f,1e6*np.array([x['value'] for x in points])/f**2,
                yerr=1e6*np.array([x['mc_error'] for x in points])/f**2,
                fmt='o-',color=color,label=f'N = {row["N"]:g}')
        ax.set_xlabel(r'$f=\alpha/\alpha_{\rm physical}$')
        ax.set_ylabel(r'$[\delta(m/100)-\delta(m)]/f^2$ [$10^{-6}$]')
        ax.set_title('Mass dependence of the residual; δ = CEEX2 − FF NLL\n'
                     'Constant: quadratic; line through zero: cubic; higher powers also possible',fontsize=10)
        ax.grid(alpha=.2);ax.legend()
        for suffix in ('png','pdf','svg'):fig.savefig(a.output/('mass-coupling-scaling.'+suffix),dpi=180)
        plt.close(fig)
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    for r in rows:
        print('N',r['N'],'quadratic+cubic',r['fits']['quadratic_cubic'],
              'mass change',r['reduced_mass_change'])


if __name__=='__main__':main()
