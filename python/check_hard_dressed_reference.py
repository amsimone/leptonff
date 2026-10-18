#!/usr/bin/env python3
"""Independent Poisson, signed hard convolution and angular dressing checks."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform

from check_ll_reference import read_run, poisson_reference
from check_matching_reference import acceptance, bilinear
from check_fiducial_reference import gauss, band_quadrature, check_close
from check_reference import simpson, source_digest


def convolution(left, right):
    count = len(left)
    result = [0.0]*count
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[(i+j) % count] += x*y
    return result


def independent_projection(run):
    metadata, weights, _ = poisson_reference(run, return_data=True)
    solution = run["solutions"]["baseline"]
    count, h, cutoff = solution["bins"], solution["spacing"], solution["soft_log_cutoff"]
    log_scale = 2*math.log(run["parent_mass"]/run["mu_frag"])
    hard = [0.0]*count
    for j in range(1, count):
        lo, hi = max(cutoff, (j-.5)*h), (j+.5)*h
        if lo < hi:
            hard[j] = simpson(lambda u: 2*math.exp(-u)*(1+math.exp(-2*u))
                /(-math.expm1(-u))*(log_scale-u+2*math.log(-math.expm1(-u))),
                lo, hi, 1e-12)
    hard[0] = 1.5-sum(hard)
    evolved = convolution(hard, weights)
    for name, values in (("initial_hard_head", hard), ("evolved_hard_head", evolved)):
        check_close(max(abs(x-y) for x, y in zip(values, solution[name])), 0,
                    2e-9, f"independent {name}")
    fractions = [math.exp(-j*h) for j in range(count)]
    quadrature_moment_residuals = []
    for power in range(4):
        quadrature = sum(node["weight"]*node["fraction"]**power
                         for node in solution["quadrature_nodes"])
        exact_lattice = sum(w*z**power for w, z in zip(weights, fractions))
        residual = quadrature-exact_lattice
        check_close(residual, 0, 2e-12, "independent compressed LL moment")
        quadrature_moment_residuals.append(residual)
    bare = sum(x*acceptance(run, z, 1) for x, z in zip(hard, fractions))
    dressed = bilinear(run, evolved, weights, fractions)
    shift = dressed-bare
    check_close(shift, solution["projected_dressing_shift"], 2e-8,
                "independent signed hard projection dressing")
    # Pure mass-window convolution reduces to hard x U x U, independently
    # of the two-dimensional prefix-sum measurement used by the driver.
    if run["pt_min"] == 0 and math.tanh(run["eta_max"]) == 1:
        combined = convolution(hard, convolution(weights, weights))
        reduced = sum(w for j, w in enumerate(combined)
                      if run["mass_min"] < run["parent_mass"]*math.exp(-j*h/2)
                      < run["mass_max"])
        check_close(reduced, dressed, 2e-8, "independent mass-window U(2t) reduction")
    metadata["quadrature_moment_residuals"] = quadrature_moment_residuals
    metadata["hard_projection_shift_residual"] = shift-solution["projected_dressing_shift"]
    return metadata


def pair_geometry(run, f1, f2, count):
    s = run["solutions"]["baseline"]
    mass, pt, eta = run["parent_mass"], run["pt_min"], run["eta_max"]
    if f1*f2 == 0 or min(f1, f2) <= 2*pt/mass:
        return 0.0
    lo = (run["mass_min"]/mass)**2/(f1*f2)
    hi = min(1-s["hard_soft_cutoff"], (run["mass_max"]/mass)**2/(f1*f2))
    if lo >= hi:
        return 0.0
    eta_bound = math.tanh(eta)

    def bound(energy):
        if energy <= pt:
            return 0.0
        return min(eta_bound, math.sqrt((1-pt/energy)*(1+pt/energy)))

    def outer(z):
        b12, b21 = min(bound(mass*z*f1/2), bound(mass*f2/2)), \
                    min(bound(mass*z*f2/2), bound(mass*f1/2))
        def angular(u):
            c = 1-u*u
            x1, x2 = (1+z-(1-z)*c)/2, (1+z+(1-z)*c)/2
            cosine = min(1.0, abs(x1*x2-2*z)/(x1*x2))
            sine = math.sqrt(z)*(1-z)*math.sqrt((1-c)*(1+c))/(x1*x2)
            a12 = band_quadrature(bound(mass*x1*f1/2), bound(mass*x2*f2/2), cosine, sine)
            a21 = band_quadrature(bound(mass*x1*f2/2), bound(mass*x2*f1/2), cosine, sine)
            density = (x1*x1+x2*x2)/((1-x1)*(1-x2))*(1-z)/2
            return 2*u*density*(a12-b12+a21-b21)
        edges = [math.sqrt(s["hard_collinear_cutoff"]), 1.0]
        for fraction in (f1, f2):
            c = ((1+z)-4*pt/(mass*fraction))/(1-z)
            if 0 < c < 1:
                u = math.sqrt(1-c)
                if edges[0] < u < 1:
                    edges.append(u)
        edges.sort()
        return sum(gauss(angular, a, b, count) for a, b in zip(edges, edges[1:]))

    edges = sorted(set([lo, hi]+[z for z in (
        2*pt/(mass*f1), 2*pt/(mass*f2), 2*pt/mass*(1/f1+1/f2)-1) if lo < z < hi]))
    result = 0.0
    for a, b in zip(edges, edges[1:]):
        def mapped(v):
            sine, cosine = math.sin(v), math.cos(v)
            return outer(a+(b-a)*sine*sine)*2*(b-a)*sine*cosine
        result += gauss(mapped, 0, math.pi/2, count)
    return result


def independent_geometry(run, count):
    nodes = run["solutions"]["baseline"]["quadrature_nodes"]
    shift = 0.0
    for i, left in enumerate(nodes):
        for j in range(i, len(nodes)):
            right = nodes[j]
            weight = left["weight"]*right["weight"]*(1 if i == j else 2)
            if i == 0 and j == 0:
                weight -= 1
            if abs(weight) < 1e-13:
                continue
            shift += weight*pair_geometry(run, left["fraction"], right["fraction"], count)
    return shift


def verify_diagnostics(run):
    a = run["alpha"]/(2*math.pi)
    for name, s in run["solutions"].items():
        for field in ("hard_number_residual", "hard_imaginary_weight_mass"):
            check_close(s[field], 0, 1e-7, f"{name} {field}")
        for field in ("quadrature_number_residual", "quadrature_first_moment_residual"):
            check_close(s[field], 0, 1e-12, f"{name} {field}")
        check_close(s["dressed_hard_coefficient"],
                    s["undressed_hard_coefficient"]+s["projected_dressing_shift"]
                    +s["geometry_dressing_shift"], 1e-11, "full hard coefficient")
        check_close(s["matched_rate"], s["ff_matched_rate"]+a*s["dressed_hard_coefficient"],
                    1e-12, "one finite hard insertion")
        if a*(s["hard_wrap_coefficient_bound"]+2*s["hard_omitted_tail_coefficient_bound"]) > 1e-10:
            raise RuntimeError("Hard projection wrap/tail bound exceeded")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--input", type=Path, help="Existing production-driver JSON, avoiding a second run")
    parser.add_argument("--coarse-only", action="store_true", help="Fast independent CTest fixture")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    executable, config = args.executable.resolve(), args.config.resolve()
    coarse = read_run(executable, config, 128, .15, 32, 4)
    result = (json.loads(args.input.read_text()) if args.input else
              coarse if args.coarse_only else read_run(executable, config))
    for field in ("library_version", "compiler", "real_digits", "model", "scheme", "coupling",
                  "parent_mass", "muon_mass", "mu_0", "mu_frag", "alpha", "pt_min", "eta_max",
                  "mass_min", "mass_max", "absolute_tolerance", "relative_tolerance"):
        if result[field] != coarse[field]:
            raise RuntimeError(f"Production JSON and executable/configuration differ: {field}")
    verify_diagnostics(result)
    poisson = independent_projection(coarse)
    low, high = independent_geometry(coarse, 32), independent_geometry(coarse, 64)
    expected = coarse["solutions"]["baseline"]["geometry_dressing_shift"]
    check_close(high, low, 2e-4, "independent angular refinement")
    check_close(high, expected, 1e-4, "independent LL-dressed angular remainder")
    mass_config = Path(__file__).resolve().parents[1]/"configs"/"vector_mass.conf"
    mass_run = read_run(executable, mass_config, 128, .15, 32, 4)
    mass_checks = independent_projection(mass_run)
    result["independent_checks"] = {
        "coarse_poisson": poisson, "mass_window_poisson": mass_checks,
        "geometry_32": low, "geometry_64": high,
        "geometry_coefficient_residual": high-expected,
        "coarse_driver_result": dict(coarse),
        "mass_window_driver_result": mass_run,
    }
    if not args.coarse_only:
        s = result["solutions"]
        shifts = {
            "ff_quadrature": s["refined_ff_quadrature"]["matched_rate"]-s["baseline"]["matched_rate"],
            "grid": s["refined_grid"]["matched_rate"]-s["refined_ff_quadrature"]["matched_rate"],
            "LL_cutoff": s["refined_cutoff"]["matched_rate"]-s["refined_grid"]["matched_rate"],
            "hard_geometry": s["refined_hard_geometry"]["matched_rate"]-s["refined_cutoff"]["matched_rate"],
        }
        for name, shift in shifts.items():
            check_close(shift, 0, 1e-5, f"production {name} refinement")
        result["refinement_rate_shifts"] = shifts
    result["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "mass_config_sha256": hashlib.sha256(mass_config.read_bytes()).hexdigest(),
        "input_json_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest() if args.input else None,
        "python_version": platform.python_version(), "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(f"Independent dressed angular residual {high-expected:.3e}")
    print(f"Dressed matched rate {result['solutions']['refined_hard_geometry']['matched_rate']:.12f}")
    print("Independent Poisson, hard projection, mass-window and angular dressing checks passed")


if __name__ == "__main__":
    main()
