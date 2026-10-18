#!/usr/bin/env python3
"""Independent Simpson and soft-slicing check of angular-inclusive hard matching."""

import argparse
import functools
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import tempfile

from check_reference import simpson, source_digest


@functools.lru_cache(maxsize=None)
def li2(x):
    # Integral definition rather than the C++ reflection/power-series algorithm.
    if x == 1:
        return math.pi**2 / 6
    if x == -1:
        return -math.pi**2 / 12
    def transformed(u):
        if u == 0:
            return 0.0
        if u == 1:
            return 2*x
        t = (1-u)*(1+u)
        logarithm = math.log((1-x)+x*u*u) if x > .9 else math.log1p(-x*t)
        return -2*u*logarithm/t
    return simpson(transformed, 0, 1, 2e-11)


def real_integral(mass, muon, lower, upper):
    threshold = 4*(muon/mass)**2
    lower = max(lower, threshold)
    if upper <= lower:
        return 0.0
    beta2 = 1-threshold
    # Logarithmic soft map combined with a square-root threshold map.
    v_lo = math.sqrt(-math.log((1-lower)/(1-threshold)))
    v_hi = math.sqrt(-math.log((1-upper)/(1-threshold)))

    def mapped(v):
        if v == 0:
            return 0.0
        photon = (1-threshold)*math.exp(-v*v)
        z = threshold-(1-threshold)*math.expm1(-v*v)
        ratio = threshold/z
        pair_beta = math.sqrt(max(0.0, 1-ratio))
        ell = 2*math.log1p(pair_beta)-math.log(ratio)
        # hep-ph/0505236v2, eqs. (3)-(5), angular primitives.
        a_term = -(3-beta2)*z/photon*pair_beta
        b_term = ((3-beta2)*(1+beta2)/photon-2*(3-beta2)+2*photon)*ell
        c_term = -2*photon*pair_beta
        return (2*a_term+b_term+c_term)*2*v*photon

    return simpson(mapped, v_lo, v_hi, 2e-9)


def soft_virtual(mass, muon, delta):
    # hep-ph/0505236v2, eq. (13), independently of the C++ inclusive eq. (19).
    r = (muon/mass)**2
    b = math.sqrt(1-4*r)
    ell = 2*math.log1p(b)-math.log(4*r)
    bracket = ((1+b*b)/(2*b)*ell-1)*math.log(delta)-1+math.log((1+b)/2)
    bracket += ((3-2*b+2*b*b)/(4*b)-3*(1-b*b)/(4*b*(3-b*b)))*ell
    bracket += (1+b*b)/(2*b)*(
        -2*li2(b)+2*li2(-b)+li2((1+b)/2)-li2((1-b)/2)+math.pi**2/2)
    return 4*b*(1+2*r)*bracket


@functools.lru_cache(maxsize=None)
def inclusive_by_slicing(mass, muon):
    def sliced(delta):
        return soft_virtual(mass, muon, delta)+real_integral(mass, muon, 0, 1-delta)

    coarse, fine = sliced(1e-4), sliced(5e-5)
    # Cancel the leading error of the soft approximation; expose its size too.
    return 2*fine-coarse, fine-coarse


def lp_coefficient(mass, muon, low, high):
    logarithm = 2*math.log(mass/muon)
    def density(z):
        if z == 0:
            # Simpson endpoint map below makes this limit zero.
            return 0.0
        return 2*(1+z*z)/(1-z)*(logarithm+math.log(z)-1)

    def interval(lower, upper):
        if upper <= lower:
            return 0.0
        # z=u^2 treats the integrable log(z) at zero independently of C++ GK.
        return simpson(lambda u: 0.0 if u == 0 else density(u*u)*2*u,
                       math.sqrt(lower), math.sqrt(upper), 2e-9)

    if low >= 1:
        return 0.0
    if high > 1:
        return 1.5-interval(0, low)
    return interval(low, high)


def check_close(actual, expected, tolerance, name):
    if not math.isfinite(actual) or abs(actual-expected) > tolerance:
        raise RuntimeError(f"{name}: actual={actual}, expected={expected}")


def verify(output):
    mass = output["parent_mass"]
    low, high = (output["mass_min"]/mass)**2, (output["mass_max"]/mass)**2
    a = output["alpha"]/(2*math.pi)
    records = [(output["muon_mass"], output["result"])]
    records += [(row["muon_mass"], row["result"]) for row in output["mass_scan"]]
    maximum_lp_residual = maximum_massive_residual = 0.0
    slicing_shifts = []
    for muon, result in records:
        lp = lp_coefficient(mass, muon, low, high)
        inclusive, shift = inclusive_by_slicing(mass, muon)
        slicing_shifts.append(shift)
        massive = 0.0 if low >= 1 else (
            inclusive-real_integral(mass, muon, 0, low) if high > 1
            else real_integral(mass, muon, low, high))
        check_close(result["matched_coefficient"], lp, 2e-7, "independent LP coefficient")
        check_close(result["massive_coefficient"], massive, 2e-6,
                    "independent soft+virtual+real coefficient")
        check_close(result["leading_power_rate"], result["born"]+a*lp,
                    3e-10, "independent LP rate")
        check_close(result["massive_rate"], result["massive_born"]+a*massive,
                    3e-9, "independent massive rate")
        r = (muon/mass)**2
        check_close(result["massive_rate"], result["leading_power_rate"],
                    100*a*r*(1+abs(math.log(r)))+10*r*r, "mass power suppression")
        maximum_lp_residual = max(maximum_lp_residual, abs(lp-result["matched_coefficient"]))
        maximum_massive_residual = max(maximum_massive_residual,
                                       abs(massive-result["massive_coefficient"]))
    power_residuals = [abs(row["result"]["massive_rate"]-row["result"]["leading_power_rate"])
                       for row in output["mass_scan"]]
    for larger, smaller in zip(power_residuals, power_residuals[1:]):
        if smaller > .6*larger+1e-10:
            raise RuntimeError("Mass scan does not approach leading power")
    rates = [row["result"]["leading_power_rate"] for row in output["scale_scan"]]
    check_close(max(rates), min(rates), 1e-10, "hard+FF fragmentation scale cancellation")
    return {
        "max_lp_coefficient_residual": maximum_lp_residual,
        "max_massive_coefficient_residual": maximum_massive_residual,
        "soft_slicing_coefficient_shifts": slicing_shifts,
        "scale_rate_span": max(rates)-min(rates),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    executable, config = args.executable.resolve(), args.config.resolve()
    base = config.read_text()
    cases = [("configured", {})]
    # Separate real-only and inclusive windows exercise both virtual branches.
    cases += [("radiative_tail", {"mass_min": 70, "mass_max": 85}),
              ("inclusive", {"mass_min": 0, "mass_max": 120}),
              ("heavier_muon", {"muon_mass": 2, "mu_0": 2})]
    reports = []
    with tempfile.TemporaryDirectory(prefix="leptonff-hard-") as directory:
        for name, changes in cases:
            lines = []
            for line in base.splitlines():
                key = line.split("=", 1)[0].strip()
                lines.append(f"{key} = {changes[key]}" if key in changes else line)
            path = Path(directory)/f"{name}.conf"
            path.write_text("\n".join(lines)+"\n")
            output = json.loads(subprocess.run([str(executable), str(path)], check=True,
                                               capture_output=True, text=True).stdout)
            checks = verify(output)
            reports.append({"case": name, "config_text": path.read_text(),
                            "output": output, "checks": checks})
            print(f"{name}: LP={output['result']['leading_power_rate']:.12f}, "
                  f"massive={output['result']['massive_rate']:.12f}, "
                  f"independent coefficient residual={checks['max_massive_coefficient_residual']:.3e}")
    report = {
        "cases": reports,
        "provenance": {
            "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
            "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
            "python_version": platform.python_version(), "platform": platform.platform(),
        },
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print("Independent hard+FF, scale, mass and soft-slicing checks passed")


if __name__ == "__main__":
    main()
