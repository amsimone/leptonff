#!/usr/bin/env python3
"""Independent causal FF series, hard projection and fiducial tensor check."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_reference import source_digest
from check_fiducial_reference import gauss_rule, projected_full_coefficient, geometry_coefficient
from check_hard_dressed_reference import pair_geometry
from check_photon_pair_reference import fiducial as photon_geometry, hard_moment
from check_pair_moment_reference import muon_hard_moment
from check_moment_reference import analytic_running
from check_xspace_reference import independent_cells, power_evolution, convolution, unit, add, scale


def hard_cells(n, span, count):
    h = span/n
    out = [[0.]*n for _ in range(2)]

    def densities(z, q, l, s):
        return (2*(1+z*z)/q*(l+2*s),
                (1+q*q)/z*(2*s+l+1)-4*q/z+q*q/(2*z))

    for cell in range(n-1):
        za, zb = math.exp(-h*cell), math.exp(-h*(cell+1))
        for x, w in gauss_rule(count):
            t = (x+1)/2
            if cell == 0:
                q = (1-zb)*t**4
                z, l, s = 1-q, math.log1p(-q), math.log(q)
                right, jac = q/(1-zb), 2*w*(1-zb)*t**3
            else:
                z = zb+(za-zb)*t
                q, l, s = 1-z, math.log(z), math.log1p(-z)
                right, jac = (za-z)/(za-zb), w*(za-zb)/2
            for a, density in zip(out, densities(z, q, l, s)):
                v = jac*z*density
                if cell:
                    a[cell] += v*(1-right)
                a[cell+1] += v*right
    cap = math.exp(-h*(n-1))
    tail = [0., 0.]
    for x, w in gauss_rule(count):
        t = (x+1)/2
        z = cap*t**4
        for i, density in enumerate(densities(z, 1-z, math.log(z), math.log1p(-z))):
            tail[i] += 2*w*cap*t**3*z*density
    for a, moment, omitted in zip(out, (muon_hard_moment(2, 0), hard_moment(2, 0)), tail):
        a[0] = moment-sum(a[1:])-omitted
    return out


def causal(a, b):
    return [sum(a[k]*b[j-k] for k in range(j+1)) for j in range(len(a))]


def acceptance(run, z, w):
    M = run["parent_mass"]
    mass = M*math.sqrt(z*w)
    ratio = 2*run["pt_min"]/(M*min(z, w))
    if ratio >= 1 or not run["mass_min"] < mass < run["mass_max"]:
        return 0.0
    return min(math.tanh(run["eta_max"]), math.sqrt(1-ratio*ratio))


def compress(weights, fractions, minimum, panels, h):
    out = [weights[0]]+[0.]*(len(fractions)-1)
    loss = math.sqrt(1-minimum)
    for j, w in enumerate(weights[1:], 1):
        z = math.exp(-h*j)
        k = min(panels-1, int(panels*math.sqrt(1-z)/loss))
        lo, hi = 1-(loss*(k+1)/panels)**2, 1-(loss*k/panels)**2
        selected = [i for i, x in enumerate(fractions) if i and lo < x < hi]
        if len(selected) != 4:
            raise RuntimeError("Missing active cubic interpolation panel")
        for i in selected:
            basis = math.prod((z-fractions[l])/(fractions[i]-fractions[l])
                              for l in selected if l != i)
            out[i] += w*basis/z
    return out


def check(run):
    if run["mu_0"] != run["muon_mass"] or run["mu_frag"] != run["parent_mass"]:
        raise ValueError("Independent fixture requires mu_0=m_mu and mu_frag=M")
    row = run["points"][0]
    n, span = row["bins"], row["log_span"]
    h, a = span/n, row["a_initial"]
    length = len(row["hard_muon"])
    L = 2*math.log(run["mu_frag"]/run["mu_0"])
    cells = independent_cells(n, span, 64)
    cells = [[[[v for v in channel[:length]] for channel in parent] for parent in A] for A in cells]
    P0, P1, B = cells
    E = power_evolution(P0, P1, a, L, 20)
    D = convolution(E, add(unit(length), B, a))
    coarse = power_evolution(P0, P1, a, L, 16)
    c1 = add(B, P0, L)
    c2 = add(scale(P1, L), convolution(P0, B), L)
    c2 = add(c2, convolution(P0, P0), L*L/2)
    c2 = add(c2, P0, L*L/3)
    arrays = {"same": D[0][0], "opposite": D[0][1], "gamma": D[2][0],
              "c1_same": c1[0][0], "c1_opposite": c1[0][1],
              "c1_gamma": c1[2][0], "c2_same": c2[0][0]}
    err = {"ff_weights": max(abs(arrays[k][j]-row["weighted_arrays"][k][j])
                              for k in ("same", "opposite", "gamma") for j in range(length)),
           "ff_expansion": max(abs(arrays[k][j]-row["weighted_arrays"][k][j])
                                for k in arrays if k.startswith("c") for j in range(length)),
           "power_degree_shift": max(abs(E[i][j][k]-coarse[i][j][k])
                                     for i in range(3) for j in range(3) for k in range(length))}
    hm, hg = [w[:length] for w in hard_cells(n, span, 64)]
    err["hard_projection"] = max(abs(v-w) for key, a1 in (("hard_muon", hm), ("hard_photon", hg))
                                  for v, w in zip(a1, row[key]))
    minimum = max(2*run["pt_min"]/run["parent_mass"], (run["mass_min"]/run["parent_mass"])**2)
    f = row["fractions"]
    q = {key: compress(w, f, minimum, row["ff_panels"], h) for key, w in arrays.items()}
    I = [1.]+[0.]*(length-1)
    qi = [1.]+[0.]*(len(f)-1)

    def born(u, v):
        return sum(x*y*math.exp(h*(i+j))*acceptance(run, math.exp(-h*i), math.exp(-h*j))
                   for i, x in enumerate(u) for j, y in enumerate(v))

    def hp(H, u, v):
        return born(causal(H, u), v)

    geom_run = dict(run, solutions={"baseline": {"hard_soft_cutoff": run["soft_cutoff"],
                                                "hard_collinear_cutoff": run["collinear_cutoff"]}})
    phot_run = dict(run, collinear_cutoff_final=run["collinear_cutoff"])

    def geometry(count):
        gm = [[0.]*len(f) for _ in f]
        gg = [[0.]*len(f) for _ in f]
        for i, x in enumerate(f):
            for j, y in enumerate(f):
                gm[i][j] = pair_geometry(geom_run, x, y, count) if j >= i else gm[j][i]
                gg[i][j] = photon_geometry(phot_run, {"photon_fraction": x, "muon_fraction": y}, count)[1]
        return gm, gg

    gm0, gg0 = geometry(96)
    gm, gg = geometry(192)

    def tensor(u, G, v):
        return sum(x*G[i][j]*y for i, x in enumerate(u) for j, y in enumerate(v))

    def evaluate(gm, gg):
        S, O, G = (arrays[key] for key in ("same", "opposite", "gamma"))
        C, CO, CG, C2 = (arrays[key] for key in ("c1_same", "c1_opposite", "c1_gamma", "c2_same"))
        qS, qO, qG, qC, qCG = (q[key] for key in ("same", "opposite", "gamma", "c1_same", "c1_gamma"))
        plus = lambda u, v: [x+y for x, y in zip(u, v)]
        out = {"ff": born(S, S)+born(O, O),
               "muon_projected": hp(hm, S, S)+hp(hm, O, O),
               "muon_geometry": tensor(qS, gm, qS)+tensor(qO, gm, qO),
               "photon_projected": 2*hp(hg, G, plus(S, O)),
               "photon_geometry": 2*tensor(qG, gg, plus(qS, qO))}
        H0 = hp(hm, I, I)+gm[0][0]
        first = 2*born(C, I)+H0
        Lh = 2*math.log(run["parent_mass"]/run["mu_0"])
        second = (2*born(C2, I)+born(C, C)+born(CO, CO)+hp(hm, C, I)+hp(hm, I, C)
                  +2*tensor(qC, gm, qi)+(2/3)*Lh*H0+2*(hp(hg, CG, I)+tensor(qCG, gg, qi)))
        ah = analytic_running(a, Lh, "two_loop")[0]
        target_run = dict(run, hard_controls_final={"soft_cutoff": run["soft_cutoff"],
                "collinear_cutoff": run["collinear_cutoff"]})
        continuum = projected_full_coefficient(run)+geometry_coefficient(target_run, 96)
        out["matched"] = (out["ff"]+ah*sum(out[k] for k in ("muon_projected", "muon_geometry",
                                                            "photon_projected", "photon_geometry"))
                          +a*(continuum-first))
        return out, first, second

    low, _, second_low = evaluate(gm0, gg0)
    high, first, second = evaluate(gm, gg)
    _, _, algebra_second = evaluate(row["geometry_muon"], row["geometry_photon"])
    err["geometry_contraction"] = max(abs(high[k]-row[k]) for k in ("muon_geometry", "photon_geometry"))
    err["geometry_refinement"] = max(abs(high[k]-low[k]) for k in ("muon_geometry", "photon_geometry"))
    err["matched_rate"] = abs(high["matched"]-row["matched"])
    err["lattice_coefficient1"] = abs(first-row["lattice_coefficient1"])
    err["coefficient2"] = abs(second-row["coefficient2"])
    err["coefficient2_refinement"] = abs(second-second_low)
    err["coefficient2_algebra"] = abs(algebra_second-row["coefficient2"])
    print(json.dumps(err, indent=2))
    limits = {"ff_weights": 2e-9, "ff_expansion": 2e-6, "power_degree_shift": 1e-12,
              "hard_projection": 2e-7, "geometry_contraction": 1e-4, "geometry_refinement": 3e-4,
              "matched_rate": 1e-7, "lattice_coefficient1": 5e-6,
              # Independent sharp-cut geometry limits C2 to about 1e-2,
              # equivalent to 1.4e-8 in the rate. Its tensor/expansion algebra
              # is checked more tightly using the separately audited matrix.
              "coefficient2": 1e-2, "coefficient2_refinement": 3e-2,
              "coefficient2_algebra": 2e-5}
    for key, value in err.items():
        if value > limits[key]:
            raise RuntimeError(f"Independent coupled fiducial check failed: {key}")
    return err


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--input", type=Path, help="Attach a completed refinement scan to the checked coarse fixture")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    exe, config = args.executable.resolve(), args.config.resolve()
    coarse = json.loads(subprocess.check_output([str(exe), str(config), "256", "2", "--dump"], text=True))
    errors = check(coarse)
    run = json.loads(args.input.read_text()) if args.input else coarse
    for key in ("version", "compiler", "alpha_msbar_at_mu0", "muon_mass", "mu_0", "mu_frag",
                "parent_mass", "pt_min", "eta_max", "mass_min", "mass_max"):
        if run[key] != coarse[key]:
            raise RuntimeError(f"Production/coarse fixture differs: {key}")
    run["independent_checks"] = errors
    if args.input:
        run["independent_coarse"] = coarse
    run["provenance"] = {"source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(exe.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "input_json_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest() if args.input else None,
        "python_version": platform.python_version(), "platform": platform.platform()}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(run, indent=2, allow_nan=False)+"\n")
    print("Independent coupled fiducial hard matching checks passed")


if __name__ == "__main__":
    main()
