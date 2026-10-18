#!/usr/bin/env python3
"""Independent z-cell quadrature and causal power-series coupled evolution.

No FFTs and no Runge--Kutta steps. This validates the weighted lattice and its
FF-only fiducial match; it does not certify a hard-matched NLL cross section.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_fiducial_reference import gauss_rule
from check_matching_reference import bilinear
from check_moment_reference import analytic_inputs
from check_nlo_kernel_reference import terms
from check_reference import independent_coefficient, source_digest


def zeros(n):
    return [[[0.]*n for _ in range(3)] for _ in range(3)]


def unit(n):
    result = zeros(n)
    for i in range(3):
        result[i][i][0] = 1.
    return result


def add(a,b,scale=1.):
    return [[[x+scale*y for x,y in zip(ar,br)] for ar,br in zip(ai,bi)] for ai,bi in zip(a,b)]


def scale(a,factor):
    return [[[factor*x for x in row] for row in parent] for parent in a]


def error(a,b):
    return max(abs(x-y) for ai,bi in zip(a,b) for ar,br in zip(ai,bi) for x,y in zip(ar,br))


def convolution(a,b):
    n = len(a[0][0])
    result = zeros(n)
    for i in range(3):
        for j in range(3):
            out = result[i][j]
            for k in range(3):
                left, right = a[i][k], b[k][j]
                for u,x in enumerate(left):
                    if x:
                        for v in range(n-u):
                            out[u+v] += x*right[v]
    return result


def densities(z,q,l,s):
    p, g, c = (1+z*z)/q, (1+q*q)/z, z*z+q*q
    p0 = [[p,0.,g],[0.,p,g],[c,c,0.]]
    r = terms(z,q,l,s)
    same = (r["photonic"]+r["nonsinglet_pair"]+r["singlet_pair"]
            +r["same_charge_interference"]-20/(9*q))
    opposite = r["singlet_pair"]+r["opposite_charge_interference"]
    p1 = [[same,opposite,r["lepton_to_photon"]],
          [opposite,same,r["lepton_to_photon"]],
          [r["photon_to_lepton"],r["photon_to_lepton"],r["photon_to_photon"]]]
    boundary = [[p*(-2*s-1),0.,g*(-2*l-1)],
                [0.,p*(-2*s-1),g*(-2*l-1)],[0.,0.,0.]]
    return p0,p1,boundary


def independent_cells(n,span,count):
    h = span/n
    arrays = [zeros(n) for _ in range(3)]
    moments = analytic_inputs(2)
    for cell in range(n-1):
        za,zb = math.exp(-h*cell),math.exp(-h*(cell+1))
        for node,w in gauss_rule(count):
            t = (node+1)/2
            if cell == 0:
                q = (1-zb)*t**4
                z, l, s = 1-q,math.log1p(-q),math.log(q)
                right, jac = q/(1-zb),2*w*(1-zb)*t**3
            else:
                z = zb+(za-zb)*t
                q, l, s = 1-z,math.log(z),math.log1p(-z)
                right, jac = (za-z)/(za-zb),w*(za-zb)/2
            values = densities(z,q,l,s)
            for a,density in zip(arrays,values):
                for i in range(3):
                    for j in range(3):
                        value = jac*z*density[i][j]
                        if cell:
                            a[i][j][cell] += value*(1-right)
                        a[i][j][cell+1] += value*right
    cap = math.exp(-h*(n-1))
    tail = [[[0.]*3 for _ in range(3)] for _ in range(3)]
    for node,w in gauss_rule(count):
        t = (node+1)/2
        z = cap*t**4
        values = densities(z,1-z,math.log(z),math.log1p(-z))
        for out,density in zip(tail,values):
            for i in range(3):
                for j in range(3):
                    out[i][j] += 2*w*cap*t**3*z*density[i][j]
    for array,moment,omitted in zip(arrays,moments,tail):
        for i in range(3):
            for j in range(3):
                array[i][j][0] = moment[i][j]-sum(array[i][j][1:])-omitted[i][j]
    return arrays


def power_evolution(p0,p1,a0,L,degree):
    n = len(p0[0][0])
    # Coefficients in t=L'/L. Scalar running is solved as a formal series,
    # then the matrix coefficients use causal (non-periodic) convolution.
    aa, aa2 = [a0], []
    series = [unit(n)]
    conv0, conv1 = [], []
    result = unit(n)
    for m in range(degree):
        aa2.append(sum(aa[k]*aa[m-k] for k in range(m+1)))
        cubic = sum(aa2[k]*aa[m-k] for k in range(m+1))
        aa.append(L/(m+1)*(2*aa2[m]/3+cubic))
        conv0.append(convolution(p0,series[m]))
        conv1.append(convolution(p1,series[m]))
        next_term = zeros(n)
        for k in range(m+1):
            next_term = add(next_term,conv0[m-k],L*aa[k]/(m+1))
            next_term = add(next_term,conv1[m-k],L*aa2[k]/(m+1))
        series.append(next_term)
        result = add(result,next_term)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable",type=Path)
    parser.add_argument("config",type=Path)
    parser.add_argument("--coarse-only",action="store_true")
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    executable, config = args.executable.resolve(),args.config.resolve()
    command = [str(executable),str(config)]+(["64"] if args.coarse_only else [])
    run = json.loads(subprocess.run(command,check=True,capture_output=True,text=True).stdout)
    if run["mu_0"] != run["muon_mass"]:
        raise ValueError("Independent boundary quadrature requires mu_0=m_mu")
    point = run["points"][0]
    n, span = point["bins"],point["log_span"]
    if point["label"] != "independent_coarse" or n != 64:
        raise RuntimeError("Unexpected coarse diagnostic contract")
    p0,p1,b1 = independent_cells(n,span,40)
    coarse = independent_cells(n,span,24)
    a0, L = run["alpha_msbar_at_mu0"]/(2*math.pi),2*math.log(run["mu_frag"]/run["mu_0"])
    e = power_evolution(p0,p1,a0,L,20)
    e_low = power_evolution(p0,p1,a0,L,16)
    d = convolution(e,add(unit(n),b1,a0))
    c1 = add(b1,p0,L)
    c2 = add(scale(p1,L),convolution(p0,b1),L)
    c2 = add(c2,convolution(p0,p0),L*L/2)
    c2 = add(c2,p0,L*L/3)
    errors = {
        "cell_residual": max(error(a,point[key]) for a,key in zip((p0,p1,b1),("kernel0","kernel1","boundary1"))),
        "cell_refinement": max(error(a,b) for a,b in zip((p0,p1,b1),coarse)),
        "causal_evolution_residual": max(error(e,point["propagator"]),error(d,point["weights"])),
        "series_refinement":error(e,e_low),
        "expansion_residual":max(error(c1,point["coefficient1"]),error(c2,point["coefficient2"])),
    }
    fractions = [math.exp(-span*k/n) for k in range(n)]
    same = [w/z for w,z in zip(d[0][0],fractions)]
    opposite = [w/z for w,z in zip(d[0][1],fractions)]
    first = [w/z for w,z in zip(c1[0][0],fractions)]
    identity = [1.]+[0.]*(n-1)
    raw = bilinear(run,same,same,fractions)+bilinear(run,opposite,opposite,fractions)
    lattice_overlap = 2*bilinear(run,first,identity,fractions)
    born, continuum = independent_coefficient(run)
    matched = raw+a0*(continuum-lattice_overlap)
    errors["fiducial_residual"] = max(abs(raw-point["resummed_ff_pairs"]),
        abs(matched-point["matched_ff_pairs"]),abs(born-point["born"]))
    errors["continuum_coefficient_residual"] = abs(continuum-point["continuum_coefficient1"])
    if (errors["cell_residual"]>2e-7 or errors["cell_refinement"]>2e-7
        or errors["causal_evolution_residual"]>2e-9 or errors["series_refinement"]>2e-10
        or errors["expansion_residual"]>2e-6 or errors["fiducial_residual"]>2e-9
        or errors["continuum_coefficient_residual"]>2e-7):
        raise RuntimeError(f"Independent coupled x-space check failed: {errors}")
    for p in run["points"]:
        if p["weighted_wrap_bound"]>1e-10 or p["weighted_tail_bound"]>1e-10:
            raise RuntimeError("Uncontrolled weighted wrap/tail")
        if max(m["maximum_residual"] for m in p["moments"][:2])>1e-9:
            raise RuntimeError("Endpoint moment recovery failed")
    points = {p["label"]:p for p in run["points"]}
    if not (points["grid_x4"]["moments"][2]["maximum_residual"]
            < points["grid_x2"]["moments"][2]["maximum_residual"]
            < points["baseline"]["moments"][2]["maximum_residual"]):
        raise RuntimeError("Higher moment refinement did not converge")
    if abs(points["span_x2"]["matched_ff_pairs"]-points["grid_x4"]["matched_ff_pairs"])>1e-9:
        raise RuntimeError("Domain extension at fixed spacing changed the fiducial action")
    run["independent_checks"] = errors
    run["provenance"] = {
        "source_sha256":source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256":hashlib.sha256(executable.read_bytes()).hexdigest(),
        "config_sha256":hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version":platform.python_version(),"platform":platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(run,indent=2,allow_nan=False)+"\n")
    print(json.dumps(errors,indent=2))
    print("Independent coupled x-space and FF fiducial checks passed")


if __name__ == "__main__":
    main()
