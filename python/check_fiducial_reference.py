#!/usr/bin/env python3
"""Independent angular quadrature and full hard+FF/LL matching check.

The independent fixed-order integration covers the supplied reference geometry:
Born inside the mass window and eta looser than the Born pT bound.
"""
import argparse
import functools
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import independent_coefficient, simpson, source_digest


@functools.lru_cache(maxsize=None)
def gauss_rule(count):
    nodes = []
    for i in range(1, count+1):
        x = math.cos(math.pi*(i-.25)/(count+.5))
        for _ in range(30):
            p, previous = 1.0, 0.0
            for k in range(1, count+1):
                p, previous = ((2*k-1)*x*p-(k-1)*previous)/k, p
            derivative = count*(x*p-previous)/(x*x-1)
            updated = x-p/derivative
            if abs(updated-x) < 2e-15:
                x = updated
                break
            x = updated
        nodes.append((x, 2/((1-x*x)*derivative*derivative)))
    return nodes


def gauss(function, lower, upper, count):
    mid, half = (lower+upper)/2, (upper-lower)/2
    return half*sum(w*function(mid+half*x) for x, w in gauss_rule(count))


def band_quadrature(b1, b2, cosine, sine):
    # Directly integrate the allowed azimuth interval; no spherical-cap formula.
    if b1 == 0 or b2 == 0:
        return 0.0
    if b2 == 1:
        return b1
    if b1 == 1:
        return b2
    if sine < 1e-14:
        return min(b1, b2)
    if math.atan2(sine, cosine) <= abs(math.acos(b1)-math.acos(b2)):
        return min(b1, b2)
    roots = [abs(b2*cosine-sign*math.sqrt(1-b2*b2)*sine) for sign in (-1, 1)]
    edges = sorted([0.0, b1]+[t for t in roots if 0 < t < b1])

    def fraction(t):
        v = math.sqrt((1-t)*(1+t))*sine
        lo = max(-1.0, min(1.0, (-b2-t*cosine)/v))
        hi = max(-1.0, min(1.0, (b2-t*cosine)/v))
        return (math.acos(lo)-math.acos(hi))/math.pi

    result = 0.0
    for lo, hi in zip(edges, edges[1:]):
        def mapped(v):
            s, c = math.sin(v), math.cos(v)
            return fraction(lo+(hi-lo)*s*s)*2*(hi-lo)*s*c
        result += gauss(mapped, 0, math.pi/2, 16)
    return result


def geometry_coefficient(output, count):
    mass, pt, eta = (output[key] for key in ("parent_mass", "pt_min", "eta_max"))
    controls = output["hard_controls_final"]
    z_lo = (output["mass_min"]/mass)**2
    z_hi = min(1-controls["soft_cutoff"], (output["mass_max"]/mass)**2)
    if z_lo >= z_hi:
        return 0.0
    u_min = math.sqrt(controls["collinear_cutoff"])
    eta_bound = math.tanh(eta)

    def bound(energy):
        if energy <= pt:
            return 0.0
        return min(eta_bound, math.sqrt((1-pt/energy)*(1+pt/energy)))

    def outer(z):
        collinear = bound(mass*z/2)

        def angular(u):
            c = 1-u*u
            x1 = (1+z-(1-z)*c)/2
            x2 = (1+z+(1-z)*c)/2
            b1, b2 = bound(mass*x1/2), bound(mass*x2/2)
            cosine = min(1.0, abs(x1*x2-2*z)/(x1*x2))
            sine = math.sqrt(z)*(1-z)*math.sqrt((1-c)*(1+c))/(x1*x2)
            accepted = band_quadrature(b1, b2, cosine, sine)
            # Energy-fraction real matrix element, unlike C++ scalar trace form.
            density = (x1*x1+x2*x2)/((1-x1)*(1-x2))*(1-z)/2
            return 4*u*density*(accepted-collinear)

        u_edges = [u_min, 1.0]
        crossing = ((1+z)-4*pt/mass)/(1-z)
        if 0 < crossing < 1:
            edge = math.sqrt(1-crossing)
            if u_min < edge < 1:
                u_edges.append(edge)
        u_edges.sort()
        return sum(gauss(angular, lo, hi, count) for lo, hi in zip(u_edges, u_edges[1:]))

    z_edges = sorted(set([z_lo, z_hi]+[
        z for z in (2*pt/mass, 4*pt/mass-1) if z_lo < z < z_hi]))
    # Smooth both threshold endpoints of each outer panel.
    total = 0.0
    for lo, hi in zip(z_edges, z_edges[1:]):
        def mapped(v):
            s, c = math.sin(v), math.cos(v)
            return outer(lo+(hi-lo)*s*s)*2*(hi-lo)*s*c
        total += gauss(mapped, 0, math.pi/2, count)
    return total


def projected_full_coefficient(output):
    born, _ = independent_coefficient(output)  # also validates this fixture
    mass = output["parent_mass"]
    ratio = 2*output["pt_min"]/mass
    lower = max(ratio, (output["mass_min"]/mass)**2)
    logarithm = 2*math.log(mass/output["muon_mass"])
    below = simpson(lambda t: 0.0 if t == 0 else (
        -born*2*(1+t**4)/(1-t*t)*(logarithm+2*math.log(t)-1)*2*t),
        0, math.sqrt(lower), 2e-9)

    def above(u):
        if u == 0 or u == 1:
            return 0.0
        s, c = math.sin(math.pi*u/2), math.cos(math.pi*u/2)
        loss = (1-lower)*c*c
        z = 1-loss
        angular = math.sqrt(max(0.0, 1-(ratio/z)**2))
        # Rationalise the difference of angular acceptances at the soft endpoint.
        quotient = -ratio**2*(1+z)/(z*z*(angular+born))
        return 2*(1+z*z)*(logarithm+math.log(z)-1)*quotient*(1-lower)*math.pi*s*c

    return 1.5*born+below+simpson(above, 0, 1, 2e-9)


def check_close(actual, expected, tolerance, name):
    if not math.isfinite(actual) or abs(actual-expected) > tolerance:
        raise RuntimeError(f"{name}: actual={actual}, expected={expected}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    executable, config = args.executable.resolve(), args.config.resolve()
    result = json.loads(subprocess.run([str(executable), str(config)],
                                      check=True, capture_output=True, text=True).stdout)
    coarse, fine = geometry_coefficient(result, 48), geometry_coefficient(result, 96)
    production = result["fixed"]["geometry_coefficient"]
    check_close(fine, coarse, 1e-5, "independent geometry quadrature refinement")
    check_close(fine, production, 5e-6, "independent fiducial geometry")
    independent_full = projected_full_coefficient(result)+fine
    check_close(result["fixed"]["full_coefficient"], independent_full,
                5e-6, "independent full hard+FF coefficient")
    a = result["alpha"]/(2*math.pi)
    for key, solution in result["solutions"].items():
        check_close(solution["matched_rate"], solution["ff_matched_rate"]
                    +a*result["fixed"]["hard_coefficient"], 1e-12,
                    f"{key}: hard insertion counted exactly once")
    rates = [row["fixed"]["rate"] for row in result["scale_scan"]]
    check_close(max(rates), min(rates), 1e-10, "fixed-order fragmentation scale cancellation")
    previous = 1.0
    for row in result["mass_scan"]:
        difference = abs(row["massive"]["rate"]-row["fixed"]["rate"])
        if difference > .6*previous or difference > 5e-6:
            raise RuntimeError("Massive fiducial reference does not approach leading power")
        previous = difference
    for field in ("hard_soft_refinement_rate_shift", "hard_collinear_refinement_rate_shift",
                  "hard_tolerance_refinement_rate_shift", "massive_refinement_rate_shift"):
        check_close(result[field], 0, 3e-8, field)
    result["independent_checks"] = {
        "geometry_48": coarse, "geometry_96": fine,
        "geometry_coefficient_residual": fine-production,
        "full_coefficient": independent_full,
        "full_coefficient_residual": independent_full-result["fixed"]["full_coefficient"],
        "scale_rate_span": max(rates)-min(rates),
    }
    result["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version": platform.python_version(), "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(f"Independent geometry coefficient {fine:.10f} (residual {fine-production:.3e})")
    print(f"Full fiducial fixed-order rate  {result['fixed']['rate']:.12f}")
    print(f"Massive reference rate         {result['massive']['rate']:.12f}")
    print(f"Additive NLO+LL (refined)       {result['solutions']['refined_cutoff']['matched_rate']:.12f}")
    print("Independent angular, hard+FF, scale, mass and additive-matching checks passed")


if __name__ == "__main__":
    main()
