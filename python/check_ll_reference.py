#!/usr/bin/env python3
"""Check LL evolution with analytic moments and a non-FFT Poisson reference."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import simpson, source_digest


def read_run(executable, config, *options):
    run = subprocess.run([str(executable), str(config), *map(str, options)],
                         check=True, capture_output=True, text=True)
    return json.loads(run.stdout)


def poisson_reference(result, return_data=False):
    solution = result["solutions"]["baseline"]
    count, spacing = solution["bins"], solution["spacing"]
    cutoff, time = solution["soft_log_cutoff"], result["time"]
    rates = [0.0] * count
    for j in range(1, count):
        lo, hi = max(cutoff, (j-.5)*spacing), (j+.5)*spacing
        if lo < hi:
            rates[j] = simpson(
                lambda u: math.exp(-u)*(1+math.exp(-2*u))/(-math.expm1(-u)),
                lo, hi, tolerance=1e-13,
            )
    rate = sum(rates)
    mean = time*rate
    weights = [0.0] * count
    weights[0] = math.exp(-mean)
    term = weights.copy()
    bound = 1.0
    for order in range(1, 40):
        # Direct periodic convolution, no Fourier transform and independently
        # integrated cell rates. Endpoint/virtual factor is exp(-time*rate).
        next_term = [0.0] * count
        for i, weight in enumerate(term):
            if weight == 0:
                continue
            for j in range(1, count):
                next_term[(i+j) % count] += time/order * weight*rates[j]
        term = next_term
        weights = [left+right for left, right in zip(weights, term)]
        ratio = mean/(order+1)
        if ratio < 1:
            bound = sum(term)*ratio/(1-ratio)
            if bound < 1e-14:
                break
    else:
        raise RuntimeError("Independent Poisson series failed to converge")

    z = [math.exp(-j*spacing) for j in range(count)]
    acceptance = 0.0
    for i in range(count):
        for j in range(count):
            pt_ratio = 2*result["pt_min"]/(result["parent_mass"]*min(z[i],z[j]))
            if pt_ratio >= 1:
                continue
            mass = result["parent_mass"]*math.sqrt(z[i]*z[j])
            if result["mass_min"] < mass < result["mass_max"]:
                angular = min(math.tanh(result["eta_max"]), math.sqrt(1-pt_ratio*pt_ratio))
                acceptance += weights[i]*weights[j]*angular
    head_residual = max(abs(left-right) for left, right in
                        zip(weights, solution["weights_head"]))
    if (head_residual > 2e-11 or abs(acceptance-solution["acceptance"]) > 2e-11
            or abs(rate-solution["resolved_jump_rate"]) > 2e-10):
        raise RuntimeError("Independent Poisson/FFT LL reference mismatch")
    metadata = {
        "bins": count, "spacing": spacing, "soft_log_cutoff": cutoff,
        "poisson_tail_bound": bound, "weights_head_max_residual": head_residual,
        "acceptance_residual": acceptance-solution["acceptance"],
    }
    return (metadata, weights, rates) if return_data else metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    executable, config = args.executable.resolve(), args.config.resolve()
    result = read_run(executable, config)
    for solution in result["solutions"].values():
        for moment in solution["moments"]:
            N = int(moment["N"])
            P0 = 1.5-2*sum(1/k for k in range(1,N+1))+1/(N*(N+1))
            exact = math.exp(result["time"]*P0)
            if abs(moment["continuum"]-exact) > 2e-12:
                raise RuntimeError("Independent analytic LL moment mismatch")
            allowed = ((N-1)*solution["mean_log_loss_discretisation_bound"]
                       + solution["wrap_probability_bound"] + 1e-10)
            if abs(moment["numerical"]-exact) > allowed:
                raise RuntimeError("Numerical LL moment exceeds discretisation bound")
        if (abs(solution["normalisation_residual"]) > 1e-10
                or solution["negative_weight_mass"] > 1e-10
                or solution["imaginary_weight_mass"] > 1e-10):
            raise RuntimeError("LL normalisation/FFT checks failed")

    baseline = result["solutions"]["baseline"]["acceptance"]
    refined_grid = result["solutions"]["refined_grid"]["acceptance"]
    refined_cut = result["solutions"]["refined_cutoff"]["acceptance"]
    if abs(refined_grid-baseline) > 1e-5 or abs(refined_cut-refined_grid) > 2e-5:
        raise RuntimeError("LL fiducial refinement is not resolved for this reference")
    coarse = read_run(executable, config, 512, .04, 32)
    result["independent_poisson_reference"] = poisson_reference(coarse)
    result["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(f"Born acceptance          {result['born_acceptance']:.12f}")
    print(f"LL O(alpha) acceptance   {result['LL_first_order_acceptance']:.12f}")
    print(f"LL resummed (refined)    {refined_cut:.12f}")
    print(f"Grid refinement shift   {refined_grid-baseline:.3e}")
    print(f"Cutoff refinement shift {refined_cut-refined_grid:.3e}")
    print("Independent Poisson, continuum moments and LL refinement checks passed")


if __name__ == "__main__":
    main()
