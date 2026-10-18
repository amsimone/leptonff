#!/usr/bin/env python3
"""Independent hard-parent tensor contraction with analytic running/Magnus FFs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import source_digest
from check_moment_reference import (analytic_inputs, analytic_running, magnus,
                                    multiply, add, scaled, identity)
from check_photon_pair_reference import hard_moment


def muon_hard_moment(N, L):
    k = int(N)-1
    H = lambda n, p=1: sum(1/j**p for j in range(1, n+1))
    ll = 1.5-2*H(N)+1/(N*(N+1))
    logz = 2*H(k, 2)+1/(k+1)**2+1/(k+2)**2-1.25
    logv = H(k)**2+H(k, 2)+H(k+1)/(k+1)+H(k+2)/(k+2)-1.75
    return 1.5+2*(L*ll+logz+2*logv)


def pair(A, B, i, j):
    # Generic observed-charge tensor contraction, without charge reduction.
    return A[i][0]*B[j][1]+A[i][1]*B[j][0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    exe, config = args.executable.resolve(), args.config.resolve()
    run = json.loads(subprocess.check_output([str(exe), str(config)], text=True))
    if run["mu_0"] != run["muon_mass"]:
        raise ValueError("Independent boundary check requires mu_0=m_mu")
    a0 = run["alpha_msbar_at_mu0"]/(2*math.pi)
    M = run["parent_mass"]
    Lh = 2*math.log(M/run["mu_0"])
    inputs = {N: analytic_inputs(N) for N in (2, 3, 5)}
    errors = {key: 0.0 for key in ("hard_coefficients", "hard_coupling", "contributions",
                                  "expansion_coefficients", "magnus_refinement")}
    for row in run["moments"]:
        N, mode = int(row["N"]), row["coupling"]
        L = 2*math.log(row["mu_frag"]/run["mu_0"])
        Lf = 2*math.log(M/row["mu_frag"])
        P0, P1, d1 = inputs[N]
        af = analytic_running(a0, L, mode)[0]
        ah = analytic_running(a0, Lh, mode)[0]
        evolution = magnus(P0, P1, a0, af, L, mode, 128)
        coarse = magnus(P0, P1, a0, af, L, mode, 64)
        D = multiply(evolution, add(identity(), d1, a0))
        hm, hg = muon_hard_moment(N, Lf), hard_moment(N, Lf)
        ff = pair(D, D, 0, 1)
        mh = ah*hm*pair(D, D, 0, 1)
        gh = ah*hg*(pair(D, D, 0, 2)+pair(D, D, 1, 2))
        values = {"ff": ff, "muon_hard": mh, "photon_hard": gh, "matched": ff+mh+gh}
        errors["contributions"] = max(errors["contributions"],
                                        *(abs(row[key]-value) for key, value in values.items()))
        errors["hard_coefficients"] = max(errors["hard_coefficients"],
            abs(hm-row["hard_muon_coefficient"]), abs(hg-row["hard_photon_coefficient"]))
        errors["hard_coupling"] = max(errors["hard_coupling"], abs(ah-row["a_hard"]))
        errors["magnus_refinement"] = max(errors["magnus_refinement"],
            *(abs(evolution[i][j]-coarse[i][j]) for i in range(3) for j in range(3)))
        c1 = add(d1, P0, L)
        c2 = add(scaled(P1, L), multiply(P0, d1), L)
        c2 = add(c2, multiply(P0, P0), L*L/2)
        b0 = 0 if mode == "fixed" else 2/3
        c2 = add(c2, P0, b0*L*L/2)
        series = [identity(), c1, c2]

        def pair_order(i, j, order):
            return sum(pair(series[k], series[order-k], i, j) for k in range(order+1))

        first = pair_order(0, 1, 1)+hm
        second = pair_order(0, 1, 2)+hm*(pair_order(0, 1, 1)+b0*Lh)
        second += hg*(pair_order(0, 2, 1)+pair_order(1, 2, 1))
        errors["expansion_coefficients"] = max(errors["expansion_coefficients"],
            abs(first-row["coefficient1"]), abs(second-row["coefficient2"]))
    print(json.dumps(errors, indent=2))
    if errors["expansion_coefficients"] > 2e-7 or max(
            value for key, value in errors.items() if key != "expansion_coefficients") > 2e-8:
        raise RuntimeError("Independent hard-matched pair moment failed")
    run["independent_checks"] = errors
    run["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version": platform.python_version(), "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(run, indent=2, allow_nan=False)+"\n")
    print("Independent hard-matched pair moment checks passed")


if __name__ == "__main__":
    main()
