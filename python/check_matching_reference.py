#!/usr/bin/env python3
"""Independent signed-boundary convolution and FF-only overlap check."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

from check_ll_reference import read_run, poisson_reference
from check_reference import independent_coefficient, simpson, source_digest


def acceptance(result, z1, z2):
    ratio = 2*result["pt_min"]/(result["parent_mass"]*min(z1,z2))
    mass = result["parent_mass"]*math.sqrt(z1*z2)
    if ratio >= 1 or not result["mass_min"] < mass < result["mass_max"]:
        return 0.0
    return min(math.tanh(result["eta_max"]),math.sqrt(1-ratio*ratio))


def bilinear(result, left, right, fractions):
    return sum(x*y*acceptance(result,fractions[i],fractions[j])
               for i,x in enumerate(left) for j,y in enumerate(right)
               if 2*result["pt_min"] < result["parent_mass"]*min(fractions[i],fractions[j]))


def check_coarse(coarse):
    metadata, weights, rates = poisson_reference(coarse,return_data=True)
    solution = coarse["solutions"]["baseline"]
    count, spacing, cutoff = solution["bins"],solution["spacing"],solution["soft_log_cutoff"]
    L0 = 2*math.log(coarse["mu_0"]/coarse["muon_mass"])
    L = 2*math.log(coarse["mu_frag"]/coarse["mu_0"])
    a = coarse["alpha"]/(2*math.pi)
    initial = [0.0]*count
    for j in range(1,count):
        lo, hi = max(cutoff,(j-.5)*spacing),(j+.5)*spacing
        if lo < hi:
            initial[j] = simpson(
                lambda u: math.exp(-u)*(1+math.exp(-2*u))/(-math.expm1(-u))
                *(L0-1-2*math.log(-math.expm1(-u))),lo,hi,tolerance=1e-13)
    initial[0] = -sum(initial)
    evolved = [0.0]*count
    for i,weight in enumerate(weights):
        for j,coefficient in enumerate(initial):
            evolved[(i+j)%count] += weight*coefficient
    head_error = max(abs(x-y) for x,y in zip(initial,solution["initial_boundary_head"]))
    evolved_error = max(abs(x-y) for x,y in zip(evolved,solution["evolved_boundary_head"]))
    if head_error > 1e-9 or evolved_error > 1e-9:
        raise RuntimeError("Independent signed-boundary cell/convolution mismatch")
    fractions = [math.exp(-j*spacing) for j in range(count)]
    correction = 2*bilinear(coarse,weights,evolved,fractions)
    rates[0] = -sum(rates)
    overlap = 2*sum((L*rate+boundary)*acceptance(coarse,z,1)
                    for rate,boundary,z in zip(rates,initial,fractions))
    _, continuum_coefficient = independent_coefficient(coarse)
    raw = solution["acceptance"]+a*correction
    matched = raw+a*(continuum_coefficient-overlap)
    if (abs(correction-solution["evolved_boundary_coefficient"]) > 1e-8
            or abs(overlap-solution["lattice_overlap_coefficient"]) > 1e-8
            or abs(matched-solution["matched_FF_acceptance"]) > 2e-10):
        raise RuntimeError("Independent FF-only overlap/matching mismatch")
    metadata.update({
        "initial_boundary_head_residual": head_error,
        "evolved_boundary_head_residual": evolved_error,
        "matched_acceptance_residual": matched-solution["matched_FF_acceptance"],
    })
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable",type=Path)
    parser.add_argument("config",type=Path)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    executable,config = args.executable.resolve(),args.config.resolve()
    result = read_run(executable,config)
    _, coefficient = independent_coefficient(result)
    for solution in result["solutions"].values():
        if abs(solution["fixed_order_coefficient"]-coefficient) > 1e-7:
            raise RuntimeError("Full Oalpha FF coefficient differs from independent reference")
        if (abs(solution["boundary_initial_number_residual"]) > 1e-8
                or abs(solution["boundary_evolved_number_residual"]) > 1e-8
                or solution["boundary_imaginary_weight_mass"] > 1e-8):
            raise RuntimeError("Boundary number/FFT diagnostics failed")
    base = result["solutions"]["baseline"]["matched_FF_acceptance"]
    fine = result["solutions"]["refined_grid"]["matched_FF_acceptance"]
    cut = result["solutions"]["refined_cutoff"]["matched_FF_acceptance"]
    if abs(fine-base) > 1e-5 or abs(cut-fine) > 1e-5:
        raise RuntimeError("FF-only matched regulator dependence is unresolved")
    result["independent_signed_reference"] = check_coarse(read_run(executable,config,512,.04,32))
    result["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),"platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(f"Fixed-order FF acceptance   {result['fixed_order_FF_acceptance']:.12f}")
    print(f"Matched FF-only (refined)   {cut:.12f}")
    print(f"Grid refinement shift      {fine-base:.3e}")
    print(f"Cutoff refinement shift    {cut-fine:.3e}")
    print("Independent signed-boundary convolution and overlap checks passed")


if __name__ == "__main__":
    main()
