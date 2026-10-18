#!/usr/bin/env python3
"""Independent continuum audit of the COMPLETE normalized two-lepton radiator."""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile
import shutil

import numpy as np
from scipy.special import spence
from scipy.integrate import quad


@lru_cache(None)
def continuum(N, nodes=256):
    """Integer moments: rational first-order integrals and separate Gauss P1.

    p(z)*(z^(N-1)-1) = -(1+z^2)*sum(z^j,j=0..N-2).
    C=d1+h1 = [p(z)*(log(z)-1)]_+; endpoint soft logs cancel.
    No leptonff code, special functions or lattice enters this calculation.
    """
    if int(N)!=N or N<1:
        raise ValueError('Integer moments N>=1 required')
    P=C=0.
    for j in range(int(N)-1):
        for k in (j+1,j+3):
            P-=1/k
            C+=1/k+1/k**2
    x,w=np.polynomial.legendre.leggauss(nodes)
    t=(x+1)/2
    edge=.5*t**4
    raw=K=0.
    for upper in (False,True):
        z,q=(1-edge,edge) if upper else (edge,1-edge)
        logz=np.log1p(-q) if upper else np.log(z)
        logs=np.log(q) if upper else np.log1p(-z)
        # scipy.spence(z) = Li2(1-z), independent of the library's Li2.
        g=(1+z*z)/q*(2*logz*logs-2*logz**2-2*spence(z))+(1+z)*logz**2/2+2*z*logz-3*z+2
        raw+=float(np.dot(w*t**3,g))
        K+=float(np.dot(w*t**3,g*np.expm1((N-1)*logz)))
    delta=3/8-math.pi**2/2+6*1.2020569031595942854
    return dict(P=P,C=C,K=K,G1=delta+raw)


def prediction(row):
    c=continuum(row['N'])
    a=row['alpha']/(2*math.pi)
    L=2*math.log(row['Q']/row['mass'])
    Lf=2*math.log(row['scale']/row['mass'])
    ell=L-Lf
    r1=2*(L*c['P']+c['C'])
    ll2=2*Lf**2*c['P']**2+4*Lf*c['P']*(c['C']+ell*c['P'])
    nll2=ll2+2*Lf*c['K']
    evolution=math.exp(2*a*Lf*c['P'])
    ll=evolution*(1+2*a*(c['C']+ell*c['P']))
    nll=ll+2*a*a*Lf*evolution*c['K']
    finite=c['C']+ell*c['P']
    return dict(order_alpha=1+a*r1,ll=ll,nll=nll,
                quadratic=nll+2*a*a*evolution*finite**2,
                r1=r1,ll2=ll2,nll2=nll2,quadratic2=nll2+2*finite**2,
                delta=2*a*a*Lf*evolution*c['K'])


def p_primitive(z):
    return -2*math.log1p(-z)-z-z*z/2


def c_primitive(z):
    l=math.log(z)
    return (-2*l*math.log1p(-z)-2*spence(1-z)-(z+z*z/2)*l+
            z+z*z/4-p_primitive(z))


@lru_cache(None)
def cumulative_coefficients(cut):
    """Direct continuum double-plus convolutions, independent of FFT/moments."""
    p=p_primitive(cut);c=c_primitive(cut)
    def cross(primitive,endpoint,finite=False):
        def f(t):
            z=cut+(1-cut)*t*t
            if z>=1:return 0.
            outer=(1+z*z)/(1-z)*(math.log(z)-1 if finite else 1)
            return outer*(endpoint-primitive(cut/z))*2*(1-cut)*t
        return (c if finite else p)*endpoint+quad(f,0,1,epsabs=2e-9,epsrel=2e-10,limit=200)[0]
    def g(z):
        l=math.log(z)
        return ((1+z*z)/(1-z)*(2*l*math.log1p(-z)-2*l*l-2*spence(z))+
                (1+z)*l*l/2+2*z*l-3*z+2)
    k=-quad(g,0,cut,epsabs=2e-10,epsrel=2e-10,limit=200)[0]
    return dict(P=p,C=c,PP=cross(p_primitive,p),PC=cross(c_primitive,c),
                CC=cross(c_primitive,c,True),K=k)


def cumulative_prediction(row):
    c=cumulative_coefficients(row['cut'])
    L=2*math.log(row['Q']/row['mass']);Lf=2*math.log(row['scale']/row['mass'])
    a=row['alpha']/(2*math.pi)
    ll2=2*Lf**2*c['PP']+4*Lf*(c['PC']+(L-Lf)*c['PP'])
    extra=2*(c['CC']+
             2*(L-Lf)*c['PC']+(L-Lf)**2*c['PP'])
    return dict(order_alpha=1-2*a*(L*c['P']+c['C']),ll2=ll2,nll2=ll2+2*Lf*c['K'],
                quadratic2=ll2+2*Lf*c['K']+extra)


def index(rows):
    return {(r['Q'],r['mass'],r['scale'],r['alpha'],r['N']):r for r in rows}


def extract_second_order(rows,mode):
    # Extrapolate the complete matched result, not the FF-only auxiliary C2.
    small=sorted(rows,key=lambda r:r['alpha'])[:3]
    a=np.array([r['alpha']/(2*math.pi) for r in small])
    y=np.array([(r[mode]-r['order_alpha'])/aa**2 for r,aa in zip(small,a)])
    fit=np.polynomial.polynomial.polyfit(a/a.max(),y,2)
    return float(fit[0])


def audit(coarse,fine):
    old=index(coarse['rows']);new=index(fine['rows'])
    if set(old)!=set(new):
        raise ValueError('Refinements change audit points')
    checks=[]
    def check(name,actual,expected,tolerance):
        checks.append(dict(name=name,actual=float(actual),expected=float(expected),
                           tolerance=float(tolerance),passed=bool(abs(actual-expected)<=tolerance)))
    expected_G1=-13/4+math.pi**2/2-2*1.2020569031595942854
    check('source-sector first moment (not imposed)',continuum(1)['G1'],expected_G1,3e-10)
    for N in (1,2,3,5,10):
        check(f'independent P1 quadrature N={N}',continuum(N,128)['K'],continuum(N)['K'],3e-9)
    groups={}
    for key,r in new.items():
        ref=prediction(r);c=old[key];label=f'Q={r["Q"]},m={r["mass"]},mu={r["scale"]},a={r["alpha"]},N={r["N"]}'
        check(label+' continuum first order',r['order_alpha'],ref['order_alpha'],3e-10)
        for mode in ('ll','nll','quadratic'):
            shift=abs(r[mode]-c[mode])
            check(label+' complete '+mode+' moment',r[mode],ref[mode],4*shift+2e-10)
        delta=r['nll']-r['ll'];previous=c['nll']-c['ll']
        check(label+' P1 insertion',delta,ref['delta'],4*abs(delta-previous)+2e-11)
        check(label+' cumulative first order',r['bin_order_alpha'],
              cumulative_prediction(r)['order_alpha'],3e-10)
        if r['N']==1:
            check(label+' normalization',r['nll'],1,3e-10)
            check(label+' quadratic normalization',r['quadratic'],1,3e-10)
        groups.setdefault((r['Q'],r['mass'],r['scale'],r['N']),[]).append(r)
    coefficients=[]
    for key,rows in groups.items():
        Q,m,mu,N=key
        oldrows=[old[(r['Q'],r['mass'],r['scale'],r['alpha'],r['N'])] for r in rows]
        cref=cumulative_prediction(rows[0])
        for mode in ('ll','nll','quadratic'):
            def extract_bins(rs):
                return extract_second_order([dict(r,order_alpha=r['bin_order_alpha'],
                        **{mode:r['bin_'+mode]}) for r in rs],mode)
            value=extract_bins(rows);prior=extract_bins(oldrows)
            check(f'{key} cumulative {mode} second coefficient',value,cref[mode+'2'],
                  4*abs(value-prior)+.003)
        if N==1:
            continue  # Already checked directly; extracting 0/alpha^2 amplifies roundoff.
        ref=prediction(rows[0])
        entry=dict(Q=Q,mass=m,scale=mu,N=N)
        for mode in ('ll','nll','quadratic'):
            value=extract_second_order(rows,mode);prior=extract_second_order(oldrows,mode)
            check(f'{key} complete {mode} second coefficient',value,ref[mode+'2'],4*abs(value-prior)+.0003)
            entry[mode]=value
            entry[mode+'_previous']=prior
            if N>1:
                small=sorted(rows,key=lambda r:r['alpha'])[:3]
                residual=[r[mode]-r['order_alpha']-(r['alpha']/(2*math.pi))**2*value for r in small]
                # Adjacent alpha values differ by two; remainder must be cubic.
                if abs(residual[0])>1e-13:
                    check(f'{key} {mode} cubic remainder',residual[1]/residual[0],8,.2)
        coefficients.append(entry)
    # Fit the mass-log coefficients at mu=Q. At this central scale there is
    # no retained alpha^2 constant; fit it too, rather than assuming it zero.
    for Q in sorted({r['Q'] for r in coefficients}):
        for N in (2,3,5,10):
            rows=[r for r in coefficients if r['Q']==Q and r['scale']==Q and r['N']==N]
            L=np.array([2*math.log(Q/r['mass']) for r in rows])
            c=continuum(N)
            for mode in ('ll','nll','quadratic'):
                fit=np.polynomial.polynomial.polyfit(L,[r[mode] for r in rows],2)
                prev=np.polynomial.polynomial.polyfit(L,[r[mode+'_previous'] for r in rows],2)
                expected=[2*c['C']**2 if mode=='quadratic' else 0,
                          4*c['P']*c['C']+(2*c['K'] if mode!='ll' else 0),2*c['P']**2]
                for j,term in enumerate(('constant','L','L2')):
                    check(f'Q={Q},N={N},{mode} mass log {term}',fit[j],expected[j],4*abs(fit[j]-prev[j])+.003)
    return dict(schema='leptonff-complete-mass-audit-report-v1',
                normalization='conditioned source-sector shape; physical sector equivalence needs external check',
                checks=checks,passed=all(c['passed'] for c in checks),
                coefficients=coefficients)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--client',type=Path,required=True)
    p.add_argument('--output',type=Path)
    p.add_argument('--grids',nargs=2,type=int,default=[65536,131072])
    a=p.parse_args()
    temporary=tempfile.TemporaryDirectory() if a.output is None else None
    if temporary:a.output=Path(temporary.name)/'audit'
    a.output.mkdir(parents=True,exist_ok=False)
    client=a.output/'client';shutil.copy2(a.client,client)
    records=[]
    for cells in a.grids:
        data=subprocess.check_output([str(client.resolve()),str(cells)],text=True)
        (a.output/f'lattice-{cells}.json').write_text(data)
        records.append(json.loads(data))
    report=audit(*records)
    report['client_sha256']=hashlib.sha256(client.read_bytes()).hexdigest()
    (a.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    failed=[c for c in report['checks'] if not c['passed']]
    print(json.dumps(dict(checks=len(report['checks']),failed=failed),indent=2))
    if temporary:temporary.cleanup()
    if failed:raise SystemExit(1)


if __name__=='__main__':main()
