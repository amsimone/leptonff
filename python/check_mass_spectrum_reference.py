#!/usr/bin/env python3
"""Check mass-bin membership, hard projections, continuum Oalpha and closure.

Independent source cells/causal FF series and literal bin assignment replace
the production cumulative/prefix-sum implementation. Oalpha bins use the
independent massive leading-power spectrum and allowed-azimuth quadrature.
This does not independently reconstruct every resummed finite-geometry bin.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import source_digest
from check_xspace_reference import independent_cells, power_evolution, convolution, add, unit
from check_coupled_decay_reference import hard_cells, causal
from check_fiducial_reference import projected_full_coefficient, geometry_coefficient
from check_moment_reference import analytic_running


def verify(run):
    row = run["points"][0]
    if run["mu_0"] != run["muon_mass"] or run["mu_frag"] != run["parent_mass"]:
        raise ValueError("Independent source fixture requires mu_0=m_mu and mu_frag=M")
    n, span = row["bins"], row["log_span"]
    h, a = span/n, row["a_initial"]
    length = len(row["same"])
    source = independent_cells(n, span, 64)
    p0, p1, b1 = [[[[v for v in channel[:length]] for channel in parent] for parent in A] for A in source]
    L = 2*math.log(run["mu_frag"]/run["mu_0"])
    E = power_evolution(p0, p1, a, L, 20)
    D = convolution(E, add(unit(length), b1, a))
    same, opposite, gamma = D[0][0], D[0][1], D[2][0]
    hm, hg = [v[:length] for v in hard_cells(n, span, 64)]
    errors = {
        "ff_weights": max(abs(v-w) for key, values in (("same", same), ("opposite", opposite), ("gamma", gamma))
                          for v, w in zip(row[key], values)),
        "hard_projection": max(abs(v-w) for key, values in (("hard_muon", hm), ("hard_photon", hg))
                               for v, w in zip(row[key], values)),
        "binned_ff": 0.0, "binned_muon_projection": 0.0, "binned_photon_projection": 0.0,
        "continuum_coefficient1": 0.0, "continuum_refinement": 0.0, "fixed_order": 0.0,
        "density_normalisation": 0.0, "closure": 0.0,
    }
    edges = run["mass_edges"]

    def literal(u, v):
        output = [0.]*(len(edges)-1)
        M = run["parent_mass"]
        for i, x in enumerate(u):
            for j, y in enumerate(v):
                f, g = math.exp(-h*i), math.exp(-h*j)
                mass = M*math.sqrt(f*g)
                ratio = 2*run["pt_min"]/(M*min(f, g))
                if ratio >= 1:
                    continue
                accepted = min(math.tanh(run["eta_max"]), math.sqrt(1-ratio*ratio))
                # Independent literal cell overlap, no cumulative/prefix sums.
                for k, (lo, hi) in enumerate(zip(edges, edges[1:])):
                    if i+j==0:
                        fraction = float(lo < mass <= hi and run["mass_min"] < mass < run["mass_max"])
                    else:
                        center=(i+j)*h
                        loss_lo,loss_hi=2*math.log(M/hi),2*math.log(M/lo)
                        fraction=max(0.,min(center+h/2,loss_hi)-max(center-h/2,loss_lo))/h
                    output[k] += x*y*accepted*fraction/(f*g)
        return output

    plus = lambda x, y: [u+v for u, v in zip(x, y)]
    ff = plus(literal(same, same), literal(opposite, opposite))
    muon = plus(literal(causal(hm, same), same), literal(causal(hm, opposite), opposite))
    photon = [2*x for x in literal(causal(hg, gamma), plus(same, opposite))]

    def cumulative(lower, count):
        if lower >= run["parent_mass"]:
            return 0.0
        config = dict(run, mass_min=lower, hard_controls_final={
            "soft_cutoff": row["soft_cutoff"], "collinear_cutoff": row["collinear_cutoff"]})
        return projected_full_coefficient(config)+geometry_coefficient(config, count)

    coarse = [cumulative(edge, 48) for edge in edges]
    fine = [cumulative(edge, 96) for edge in edges]
    born = math.sqrt(1-(2*run["pt_min"]/run["parent_mass"])**2)
    for k, binrow in enumerate(row["spectrum"]):
        result = binrow["integral"]
        for key, values, target in (("binned_ff", ff, "ff"),
                                    ("binned_muon_projection", muon, "muon_projected"),
                                    ("binned_photon_projection", photon, "photon_projected")):
            errors[key] = max(errors[key], abs(values[k]-result[target]))
        c1 = fine[k]-fine[k+1]
        c1low = coarse[k]-coarse[k+1]
        errors["continuum_coefficient1"] = max(errors["continuum_coefficient1"], abs(c1-result["continuum_coefficient1"]))
        errors["continuum_refinement"] = max(errors["continuum_refinement"], abs(c1-c1low))
        born_bin = born if binrow["lower"] < run["parent_mass"] <= binrow["upper"] else 0.0
        errors["fixed_order"] = max(errors["fixed_order"], abs(born_bin+a*c1-result["fixed_order"]))
        errors["density_normalisation"] = max(errors["density_normalisation"],
            abs(binrow["matched_density_per_GeV"]*(binrow["upper"]-binrow["lower"])-result["matched"]))

    for point in run["points"]:
        for key, value in point["total"].items():
            summed = sum(binrow["integral"][key] for binrow in point["spectrum"])
            errors["closure"] = max(errors["closure"], abs(summed-value))
    limits = {"ff_weights": 2e-9, "hard_projection": 2e-7,
              "binned_ff": 2e-9, "binned_muon_projection": 2e-7, "binned_photon_projection": 2e-9,
              "continuum_coefficient1": 1e-5, "continuum_refinement": 2e-5,
              "fixed_order": 2e-8, "density_normalisation": 1e-14, "closure": 1e-10}
    print(json.dumps(errors, indent=2))
    for key, value in errors.items():
        if value > limits[key]:
            raise RuntimeError(f"Independent mass-spectrum check failed: {key}")
    return errors


def scan_checks(run):
    for point in run["points"]:
        for key, value in point["total"].items():
            if abs(sum(b["integral"][key] for b in point["spectrum"])-value) > 1e-10:
                raise RuntimeError(f"Bin partition failed: {point['label']} {key}")
        if point["weighted_wrap_bound"] > 1e-10 or point["weighted_tail_bound"] > 1e-10:
            raise RuntimeError("Spectrum FF domain bounds failed")
        if point["compression_moment_residual"] > 1e-12:
            raise RuntimeError("Spectrum FF compression moment check failed")
    points = {p["label"]: p for p in run["points"]}
    if "hard_slivers" not in points:
        return {}
    shifts = {}
    for label, left, right in (("ff_panels", "ff_panels_x2", "baseline"),
                               ("grid_x2", "grid_x2", "ff_panels_x2"),
                               ("grid_x4", "grid_x4", "grid_x2"),
                               ("span", "span_x2", "grid_x4"),
                               ("hard_slivers", "hard_slivers", "span_x2")):
        delta = [a["integral"]["matched"]-b["integral"]["matched"]
                 for a, b in zip(points[left]["spectrum"], points[right]["spectrum"])]
        shifts[label] = {"bins": delta, "maximum_bin_shift": max(abs(v) for v in delta),
                         "total_shift": points[left]["total"]["matched"]-points[right]["total"]["matched"]}
    for label in ("scale_half", "scale_double", "start_half", "start_double"):
        point = points[label]
        if abs(point["a_hard"]-points["hard_slivers"]["a_hard"]) > 1e-13:
            raise RuntimeError("Correlated variation changed hard-scale coupling")
        for current, central in zip(point["spectrum"], points["hard_slivers"]["spectrum"]):
            if abs(current["integral"]["continuum_coefficient1"]
                   -central["integral"]["continuum_coefficient1"]) > 5e-6:
                raise RuntimeError("Continuum first coefficient changed under scale variation")
        if label.startswith("start_"):
            a0=run["alpha_msbar_at_mu0"]/(2*math.pi)
            L0=2*math.log(point["mu_0"]/run["mu_0"])
            expected=analytic_running(a0,L0,"two_loop")[0]
            if abs(expected-point["a_initial"]) > 1e-14:
                raise RuntimeError("Starting-scale MSbar input was not transported consistently")
    return shifts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("bin_config", type=Path)
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    exe, config, bins = (v.resolve() for v in (args.executable, args.config, args.bin_config))
    coarse = json.loads(subprocess.check_output([str(exe), str(config), str(bins), "256", "1", "--dump"], text=True))
    errors = verify(coarse)
    run = json.loads(args.input.read_text()) if args.input else coarse
    for key in ("version", "compiler", "parent_mass", "muon_mass", "mu_0", "mu_frag", "alpha_msbar_at_mu0",
                "pt_min", "eta_max", "mass_min", "mass_max", "mass_edges"):
        if run[key] != coarse[key]:
            raise RuntimeError(f"Production/coarse inputs differ: {key}")
    run["independent_checks"] = errors
    run["refinement_rate_shifts"] = scan_checks(run)
    points={p["label"]:p for p in run["points"]}
    run["variation_rate_shifts"]={}
    if "hard_slivers" in points:
        central=points["hard_slivers"]
        for label in ("scale_half","scale_double","start_half","start_double"):
            point=points[label]
            run["variation_rate_shifts"][label]={
                "bins":[b["integral"]["matched"]-c["integral"]["matched"]
                        for b,c in zip(point["spectrum"],central["spectrum"])],
                "total_shift":point["total"]["matched"]-central["total"]["matched"]}
    if args.input:
        run["independent_coarse"] = coarse
    run["provenance"] = {"source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "bin_config_sha256": hashlib.sha256(bins.read_bytes()).hexdigest(),
        "input_json_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest() if args.input else None,
        "python_version": platform.python_version(), "platform": platform.platform()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(run, indent=2, allow_nan=False)+"\n")
    print("Independent bin membership, hard projections, Oalpha and correlated closure checks passed")


if __name__ == "__main__":
    main()
