#!/usr/bin/env python3
"""Independent time-like NLO QED matrix check; Python standard library only.

Uses transformed fixed Gauss quadrature, integral-defined dilogarithms and
analytic moments. The C++ calculation uses adaptive Gauss--Kronrod and a
dilogarithm series. This checks kernels, not an NLL evolution or pair selection.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_fiducial_reference import gauss_rule
from check_reference import source_digest

ZETA2 = math.pi**2/6
ZETA3 = 1.2020569031595942854


def dilog(x):
    if x == 1:
        return ZETA2
    if x == 0:
        return 0.
    if x > .5:
        return ZETA2-math.log(x)*math.log1p(-x)-dilog(1-x)
    # Li2(x) = -integral_0^1 dt log(1-x*t)/t; negative arguments use
    # the integral directly, unlike the C++ negative-argument reflection.
    return -.5*sum(w*math.log1p(-x*(u+1)/2)/((u+1)/2)
                   for u, w in gauss_rule(32))


def photon_terms(z, l, s):
    # Independent derivation: SL QED, de Florian et al. 1606.02887v3,
    # eqs. (55),(62),(68), continued by Stratmann--Vogelsang
    # hep-ph/9612250v1 eqs. (34)--(39). These are NOT SL kernels in the solver.
    x, lx, sx = 1/z, -l, s-l
    p = z*z+(1-z)**2
    g = (1+(1-z)**2)/z
    px, gx = x*x+(1-x)**2, (1+(1-x)**2)/x
    # Crossing ln^2(1-x) gives ln^2(x-1)-pi^2, including in (sx-lx)^2.
    sl_lgamma_crossed = .5*(4-9*x-(1-4*x)*lx-(1-2*x)*lx*lx+4*sx
        +px*(2*((sx-lx)**2-math.pi**2)-4*(sx-lx)-4*ZETA2+10))
    correction_lgamma = (2*g*(4*dilog(1-z)+2*s*l+s)
        +(-4*z*(3-z)*l-20+28*z-5*z*z)/(2*z))
    lepton_photon = z*sl_lgamma_crossed+correction_lgamma

    sl_gammal_crossed = (-(3*sx+sx*sx-math.pi**2)*gx
        +(2+3.5*x)*lx-(1-x/2)*lx*lx-2*x*sx-3.5*x-2.5
        -(4*x/3+gx*(20/9+4*sx/3)))
    correction_gammal = (4*p*(1-2*l)/3
        +2*p*(-4*dilog(1-z)-2*s*l+s)
        +(-4*(2-z)*(1-2*z)*l-5+28*z-20*z*z)/2)
    photon_lepton = z*sl_gammal_crossed+correction_gammal

    sl_gammagamma_crossed = (-16+8*x+20*x*x/3+4/(3*x)
        -(6+10*x)*lx-2*(1+x)*lx*lx)
    correction_gammagamma = 4/(9*z)*((1-z)*(38+47*z+38*z*z)
        +3*(1+z)*(4+11*z+4*z*z)*l)
    photon_photon = -z*sl_gammagamma_crossed+correction_gammagamma
    return {"lepton_to_photon": lepton_photon, "photon_to_lepton": photon_lepton,
            "photon_to_photon": photon_photon}


def terms(z, q, l, s):
    p = (1+z*z)/q
    result = {
        "photonic": p*(2*l*s-2*l*l-2*dilog(q))
                    +(1+z)*l*l/2+2*z*l-3*z+2,
        "nonsinglet_pair": -2*p*l/3-2/9+22*z/9,
        "singlet_pair": (1+z)*l*l+(-5-9*z-8*z*z/3)*l
                        -8-20/(9*z)+4*z+56*z*z/9,
        "same_charge_interference": p*(2*dilog(q)+1.5*l)
                                    -3.5*(1+z)*l-7+8*z,
        "opposite_charge_interference": 2*(1+z*z)/(1+z)
            *(-2*dilog(-z)-2*l*math.log1p(z)+l*l/2-ZETA2)
            +2*(1+z)*l+4*q,
    }
    result.update(photon_terms(z,l,s))
    return result


def integrate_rows(count, lower, N):
    result = dict.fromkeys(("photonic", "nonsinglet_pair", "singlet_pair",
                            "same_charge_interference", "opposite_charge_interference",
                            "lepton_to_photon", "photon_to_lepton", "photon_to_photon"), 0.)
    # Fourth-power maps resolve integrable log^2(z) and log(1-z) endpoints.
    # Keep q=1-z separately, so very small q never rounds to zero.
    halves = ((lower, (1+lower)/2), ((1+lower)/2, 1.))
    for index, (lo, hi) in enumerate(halves):
        for u, w in gauss_rule(count):
            t = (u+1)/2
            if index == 0:
                z = lo+(hi-lo)*t**4
                q = 1-z
                l, s = math.log(z), math.log1p(-z)
            else:
                q = (hi-lo)*t**4
                z = 1-q
                l, s = math.log1p(-q), math.log(q)
            factor = 2*w*(hi-lo)*t**3*math.exp((N-1)*l)
            values = terms(z,q,l,s)
            for name in result:
                # The singlet N1 is deliberately never integrated.
                if N != 1 or name not in ("singlet_pair","photon_to_photon"):
                    result[name] += factor*values[name]

    deltas = {"photonic": 3/8-3*ZETA2+6*ZETA3,
              "nonsinglet_pair": -4*ZETA2/3-1/6, "photon_to_photon": -1.}
    result = {name: value+deltas.get(name,0.) for name, value in result.items()}
    if lower == 0:
        result["nonsinglet_pair"] += 20/9*harmonic(N-1)
    else:
        if N != 3:
            raise ValueError("Finite-bound reference is for z^2 only")
        result["nonsinglet_pair"] += -20/9*(-1.5+lower+lower*lower/2
                                                       +math.log1p(-lower))
    # Cancel the singlet before constructing the finite N1 valence action.
    result["valence"] = (result["photonic"]+result["nonsinglet_pair"]
                         +result["same_charge_interference"]
                         -result["opposite_charge_interference"])
    result["same_charge"] = sum(result[k] for k in
        ("photonic","nonsinglet_pair","singlet_pair","same_charge_interference"))
    result["opposite_charge"] = result["singlet_pair"]+result["opposite_charge_interference"]
    if lower == 0 and N == 1:
        for name in ("singlet_pair","same_charge","opposite_charge","photon_to_photon"):
            result[name] = None
    return result


def harmonic(n, power=1):
    return sum(1/j**power for j in range(1,n+1))


def diagonal_nonsinglet_moment(N):
    # Independent log-only density, Blumlein--Kawamura eq. (30), with
    # its delta and plus coefficients reconstructed from Arbuzov--Melnikov.
    h, h2, h3 = (harmonic(N-1,k) for k in (1,2,3))
    def j(k):
        return harmonic(k)/k**2-(ZETA2-harmonic(k,2))/k
    a = 2*(h*(ZETA2-h2)+ZETA3-h3)-j(N)-j(N+1)
    b = 4*ZETA3-2*h3-2*harmonic(N+1,3)
    c = -2*ZETA2+h2+harmonic(N+1,2)
    photon_plus_interference = (3/8-3*ZETA2+6*ZETA3+2*a-2*b+1.5*c
        +1/N**3+1/(N+1)**3+7/(2*N**2)+3/(2*(N+1)**2)-5/N+5/(N+1))
    return photon_plus_interference+ns_moment(N)


def ns_moment(N):
    return (20*harmonic(N-1)/9-2*(harmonic(N-1,2)+harmonic(N+1,2))/3
            -1/6-2/(9*N)+22/(9*(N+1)))


def singlet_moment(N):
    return (2/N**3+2/(N+1)**3+5/N**2+9/(N+1)**2+8/(3*(N+2)**2)
            -8/N-20/(9*(N-1))+4/(N+1)+56/(9*(N+2)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve()
    result = json.loads(subprocess.run([str(executable)], check=True,
                        capture_output=True,text=True).stdout)
    indices = result["moment_indices"]
    if indices != [1,2,3,5] or not result["photon_channels_available"]:
        raise RuntimeError("Unexpected kernel diagnostic contract")
    if result["maximum_moment_refinement_shift"] > 2e-8:
        raise RuntimeError("C++ kernel moments are not converged")
    references = {N: integrate_rows(80,0.,N) for N in indices}
    coarse = {N: integrate_rows(48,0.,N) for N in indices}
    finite = integrate_rows(80,.37,3)
    finite_coarse = integrate_rows(48,.37,3)
    k = 13/4-3*ZETA2+2*ZETA3
    analytic = {"photonic": -k, "nonsinglet_pair": 0.,
                "same_charge_interference": 2*k,
                "opposite_charge_interference": k, "valence": 0.}
    residuals, refinements, analytic_residuals = [], [], []
    expected_names = set(finite)
    if {row["name"] for row in result["kernels"]} != expected_names:
        raise RuntimeError("Missing or unexpected kernel rows")
    for row in result["kernels"]:
        name = row["name"]
        if len(row["moments"]) != len(indices):
            raise RuntimeError("Unexpected moment count")
        for N, actual in zip(indices,row["moments"]):
            if N == 1 and name in ("singlet_pair","same_charge","opposite_charge","photon_to_photon"):
                if actual is not None:
                    raise RuntimeError("Divergent N1 must be null")
                continue
            residuals.append(actual-references[N][name])
            refinements.append(references[N][name]-coarse[N][name])
            if N == 1 and name in analytic:
                analytic_residuals.append(actual-analytic[name])
            if name == "nonsinglet_pair":
                analytic_residuals.append(actual-ns_moment(N))
            if name == "singlet_pair":
                analytic_residuals.append(actual-singlet_moment(N))
        residuals.append(row["z2_above_0p37"]-finite[name])
        refinements.append(finite[name]-finite_coarse[name])
    for N in indices:
        computed = sum(references[N][name] for name in
                       ("photonic","nonsinglet_pair","same_charge_interference"))
        analytic_residuals.append(computed-diagonal_nonsinglet_moment(N))
    for N in (2,3,5):
        expected = (-1-4/N+12/(N+1)-164/(9*(N+2))+92/(9*(N-1))
            -10/N**2-14/(N+1)**2-16/(3*(N+2)**2)-16/(3*(N-1)**2)
            +4/N**3+4/(N+1)**3)
        analytic_residuals.append(references[N]["photon_to_photon"]-expected)
    r = references[2]
    analytic_residuals.extend((r["same_charge"]+r["opposite_charge"]+r["lepton_to_photon"],
                              2*r["photon_to_lepton"]+r["photon_to_photon"]))
    if result["maximum_parent_momentum_residual"] > 2e-8:
        raise RuntimeError("NLO per-parent momentum sum failed")
    diagnostics = {
        "maximum_independent_residual": max(map(abs,residuals)),
        "maximum_independent_refinement": max(map(abs,refinements)),
        "maximum_analytic_residual": max(map(abs,analytic_residuals)),
    }
    if max(diagnostics.values()) > 2e-8:
        raise RuntimeError(f"Independent kernel check failed: {diagnostics}")
    result.update(diagnostics)
    result["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(diagnostics,indent=2))
    print("Independent complete NLO QED matrix and valence checks passed")


if __name__ == "__main__":
    main()
