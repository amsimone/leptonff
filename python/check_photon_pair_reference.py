#!/usr/bin/env python3
"""Check the photon-parent coefficient with energy variables/azimuth quadrature.

The massive check integrates the existing paper's current directly in log(t),
without the C++ analytic primitive. Hard moments use harmonic sums. Fiducial
angles use allowed-azimuth integration rather than spherical-cap intersections.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import source_digest
from check_fiducial_reference import gauss, band_quadrature


def massive_current_integral(M, m, y, count):
    r = (m/M)**2
    v = 1-y
    root = math.sqrt(v*v-4*r)
    lo = 2*r*y/(v-2*r+root)
    hi = y*(v-2*r+root)/(2*(y+r))

    def current(logt):
        t = math.exp(logt)
        z = 1-y-t
        photon = y+t
        beta = math.sqrt(1-4*r/z)
        # hep-ph/0505236v2 eq.(3), with its energy-coordinate Jacobian.
        # Use invariant expressions for 1-beta^2*c^2 near the massive pole.
        w = (y-t)/photon
        denominator = 4*y*t/photon**2
        beta2 = 1-4*r
        A = -(3-beta2)*4*r/2
        B = (3-beta2)*(1+beta2)/photon-2*(3-beta2)+2*photon
        density = beta*(2*A*(1+w*w)/(photon*denominator**2)
                        +B/denominator-photon)
        return density*2*t/(photon*beta)

    return gauss(current, math.log(lo), math.log(hi), count)


def hard_moment(N, L):
    k = int(N)-1
    harmonic = lambda n: sum(1/j for j in range(1, n+1))
    p = 2/k-2/(k+1)+1/(k+2)
    logy = -2/k**2+2/(k+1)**2-1/(k+2)**2
    logv = -2*harmonic(k)/k+2*harmonic(k+1)/(k+1)-harmonic(k+2)/(k+2)
    rational = -3.5/k+3/(k+1)+.5/(k+2)
    return (L+1)*p+logy+2*logv+rational


def fiducial(output, row, count):
    M, pt, eta = (output[k] for k in ("parent_mass", "pt_min", "eta_max"))
    fg, fm = row["photon_fraction"], row["muon_fraction"]
    lower = (output["mass_min"]/M)**2/(fg*fm)
    upper = min(1.0, (output["mass_max"]/M)**2/(fg*fm))
    if lower >= upper or 2*pt/M >= min(fg, fm):
        return 0.0, 0.0
    bgmax = math.tanh(eta)
    p, q = 2*pt/(M*fg), 2*pt/(M*fm)
    L = 2*math.log(M/output["mu_frag"])

    def bound(x, threshold):
        return min(bgmax, math.sqrt(max(0, 1-(threshold/x)**2))) if x > threshold else 0.0

    def B(y):
        return min(bound(y, p), bound(1, q))

    def hard(y):
        # Evaluate the mass-regulated current at a small auxiliary mass and
        # subtract its universal FF. This independently recovers the finite
        # MSbar coefficient; no copy of the C++ hard-spectrum formula.
        # Keep the auxiliary mass small compared with the endpoint distance;
        # a fixed mass would remove a shrinking slice near y=1.
        m = M*1e-7*min(y, 1-y)
        current = massive_current_integral(M, m, y, 32)
        g = (1+(1-y)**2)/y
        ff = g*(2*math.log(output["mu_frag"]/m)-2*math.log(y)-1)
        return (current-ff)*B(y)

    def geometry(y):
        span = 1-y
        reference = B(y)
        umin = math.asin(math.sqrt(output["collinear_cutoff_final"]))
        # t=(1-y)*sin(u)^2 smooths the opposite end of the triangle too.
        edges = [umin, math.pi/2]
        for t in (p-y, 1-q):
            if 0 < t < span:
                u = math.asin(math.sqrt(t/span))
                if umin < u < math.pi/2:
                    edges.append(u)
        edges.sort()

        def integrand(u):
            s, c = math.sin(u), math.cos(u)
            t = span*s*s
            xg, xm = y+t, 1-t
            bg, bm = bound(xg, p), bound(xm, q)
            cosine = min(1.0, abs(1-2*y/(xg*xm)))
            sine = 2*math.sqrt(max(0, y*t*(1-y-t)))/(xg*xm)
            accepted = band_quadrature(bg, bm, cosine, sine)
            xother = 1-y
            matrix = (xother*xother+xm*xm)/(y*t)
            return matrix*(accepted-reference)*2*span*s*c

        return sum(gauss(integrand, lo, hi, count) for lo, hi in zip(edges, edges[1:]))

    yedges = sorted(set([lower, upper]+[v for v in (
        p, p*math.cosh(eta), fm/fg, p+q-1) if lower < v < upper]))

    def integrate_y(function):
        def panel(lo, hi):
            def mapped(u):
                s, c = math.sin(u), math.cos(u)
                return function(lo+(hi-lo)*s*s)*2*(hi-lo)*s*c
            return gauss(mapped, 0, math.pi/2, count)
        return sum(panel(lo, hi) for lo, hi in zip(yedges, yedges[1:]))

    return integrate_y(hard), integrate_y(geometry)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    exe, config = args.executable.resolve(), args.config.resolve()
    run = json.loads(subprocess.check_output([str(exe), str(config)], text=True))
    errors = {"massive_current": 0.0, "massive_quadrature_shift": 0.0,
              "hard_moments": 0.0, "fiducial_projected": 0.0,
              "fiducial_geometry": 0.0, "fiducial_quadrature_shift": 0.0}
    for row in run["density_scan"]:
        value = massive_current_integral(run["parent_mass"], row["muon_mass"], row["y"], 64)
        refined = massive_current_integral(run["parent_mass"], row["muon_mass"], row["y"], 96)
        errors["massive_current"] = max(errors["massive_current"], abs(refined-row["massive"]))
        errors["massive_quadrature_shift"] = max(errors["massive_quadrature_shift"], abs(refined-value))
    L = 2*math.log(run["parent_mass"]/run["mu_frag"])
    for row in run["moments"]:
        errors["hard_moments"] = max(errors["hard_moments"], abs(row["hard"]-hard_moment(row["N"], L)))
    for row in run["fiducial"]:
        coarse = fiducial(run, row, 48)
        fine = fiducial(run, row, 96)
        row["independent"] = {"projected": fine[0], "geometry": fine[1],
                              "total": sum(fine)}
        errors["fiducial_projected"] = max(errors["fiducial_projected"], abs(fine[0]-row["final"]["projected"]))
        errors["fiducial_geometry"] = max(errors["fiducial_geometry"], abs(fine[1]-row["final"]["geometry"]))
        errors["fiducial_quadrature_shift"] = max(errors["fiducial_quadrature_shift"],
                                                   *(abs(a-b) for a, b in zip(coarse, fine)))
    limits = {"massive_current": 1e-9, "massive_quadrature_shift": 1e-9,
              "hard_moments": 1e-8, "fiducial_projected": 3e-6,
              "fiducial_geometry": 3e-6, "fiducial_quadrature_shift": 1e-5}
    print(json.dumps(errors, indent=2))
    for key, limit in limits.items():
        if errors[key] > limit:
            raise RuntimeError(f"Independent photon-pair check failed: {key}")
    run["independent_check"] = errors
    run["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version": platform.python_version(), "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(run, indent=2, allow_nan=False)+"\n")
    print("Independent photon-parent hard checks passed")


if __name__ == "__main__":
    main()
