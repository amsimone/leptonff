#!/usr/bin/env python3
"""Independent coupled-moment reference using analytic running and Magnus steps.

Standard library only. LO/boundary moments are analytic; NLO moments use the
independent source/continuation quadrature in check_nlo_kernel_reference.py.
The C++ reference instead evolves the coupling and matrix jointly with RK4.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from check_nlo_kernel_reference import harmonic, integrate_rows
from check_reference import source_digest


def identity():
    return [[float(i == j) for j in range(3)] for i in range(3)]


def add(a, b, scale=1.):
    return [[a[i][j]+scale*b[i][j] for j in range(3)] for i in range(3)]


def scaled(a, factor):
    return [[x*factor for x in row] for row in a]


def multiply(a, b):
    return [[sum(a[i][k]*b[k][j] for k in range(3)) for j in range(3)] for i in range(3)]


def norm(a):
    return max(abs(x) for row in a for x in row)


def difference(a, b):
    return norm(add(a,b,-1))


def exponential(a):
    halvings = max(0,math.ceil(math.log2(max(norm(a)/.25,1.))))
    a = scaled(a,2.**(-halvings))
    result, term = identity(), identity()
    for k in range(1,100):
        term = scaled(multiply(a,term),1/k)
        result = add(result,term)
        if norm(term) < 2e-17:
            break
    else:
        raise RuntimeError("Independent matrix exponential did not converge")
    for _ in range(halvings):
        result = multiply(result,result)
    return result


def analytic_inputs(N):
    ll = 1.5-2*harmonic(N)+1/(N*(N+1))
    lg = 2/(N-1)-2/N+1/(N+1)
    gl = 1/N-2/(N+1)+2/(N+2)
    p0 = [[ll,0.,lg],[0.,ll,lg],[gl,gl,-2/3]]
    logarithmic_plus = (harmonic(N-1)**2+harmonic(N-1,2)
                       +harmonic(N)/N+harmonic(N+1)/(N+1)-7/4)
    dll = -2*logarithmic_plus-ll
    dlg = 4/(N-1)**2-4/N**2+2/(N+1)**2-lg
    d1 = [[dll,0.,dlg],[0.,dll,dlg],[0.,0.,0.]]
    r = integrate_rows(80,0.,N)
    p1 = [[r["same_charge"],r["opposite_charge"],r["lepton_to_photon"]],
          [r["opposite_charge"],r["same_charge"],r["lepton_to_photon"]],
          [r["photon_to_lepton"],r["photon_to_lepton"],r["photon_to_photon"]]]
    return p0,p1,d1


def analytic_running(a0, L, mode):
    if mode == "fixed":
        return a0,a0*L,a0*a0*L
    b0, b1 = 2/3, float(mode == "two_loop")
    if mode == "one_loop":
        final = a0/(1-b0*a0*L)
        return final,-math.log1p(-b0*a0*L)/b0,(final-a0)/b0

    def elapsed(ratio):
        change = a0*(ratio-1)
        return ((ratio-1)/(b0*a0*ratio)
                +b1/(b0*b0)*(math.log1p(b1*change/(b0+b1*a0))-math.log(ratio)))
    lo, hi = (1.,1/(2*math.pi*a0)) if L >= 0 else (1e-12,1.)
    if elapsed(lo) > L or elapsed(hi) < L:
        raise RuntimeError("Requested coupling outside independent running bracket")
    for _ in range(100):
        mid = (lo+hi)/2
        if elapsed(mid) < L:
            lo = mid
        else:
            hi = mid
    final = a0*(lo+hi)/2
    beta_log = math.log1p((final-a0)/(b0+a0))
    return final,(math.log(final/a0)-beta_log)/b0,beta_log


def magnus(p0,p1,a0,final,L,mode,count):
    if mode == "fixed":
        return exponential(add(scaled(p0,a0*L),scaled(p1,a0*a0*L)))
    h = math.log(final/a0)/count
    b0, b1 = 2/3, float(mode == "two_loop")
    result = identity()
    def generator(s):
        a = a0*math.exp(s)
        return scaled(add(p0,p1,a),1/(b0+b1*a))
    for k in range(count):
        mid = (k+.5)*h
        early = generator(mid-h*math.sqrt(3)/6)
        late = generator(mid+h*math.sqrt(3)/6)
        commutator = add(multiply(early,late),multiply(late,early),-1)
        omega = add(scaled(add(early,late),h/2),commutator,-math.sqrt(3)*h*h/12)
        result = multiply(exponential(omega),result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable",type=Path)
    parser.add_argument("--output",type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve()
    result = json.loads(subprocess.run([str(executable)],check=True,
                        capture_output=True,text=True).stdout)
    if result["mu_0"] != result["muon_mass"]:
        raise ValueError("Analytic boundary reference requires mu_0=m_mu")
    expected = {(N,mode) for N in (2,3,5) for mode in ("fixed","one_loop","two_loop")}
    if {(r["N"],r["coupling"]) for r in result["moments"]} != expected or len(result["moments"]) != 9:
        raise RuntimeError("Unexpected moment diagnostic rows")
    a0 = result["alpha_msbar_at_mu0"]/(2*math.pi)
    L = 2*math.log(result["mu_frag"]/result["mu_0"])
    inputs = {N: analytic_inputs(N) for N in (2,3,5)}
    errors = dict.fromkeys(("coupling_residual","kernel_boundary_residual",
                           "evolution_residual","coefficient_residual",
                           "magnus_refinement","momentum_residual"),0.)
    for row in result["moments"]:
        N, mode = row["N"],row["coupling"]
        p0,p1,d1 = inputs[N]
        for key,reference in (("kernel0",p0),("kernel1",p1),("boundary1",d1)):
            errors["kernel_boundary_residual"] = max(errors["kernel_boundary_residual"],
                                                        difference(row[key],reference))
        final, integral_a, integral_a2 = analytic_running(a0,L,mode)
        for key,reference in (("a_final",final),("integral_a",integral_a),("integral_a2",integral_a2)):
            errors["coupling_residual"] = max(errors["coupling_residual"],abs(row[key]-reference))
        e = magnus(p0,p1,a0,final,L,mode,128)
        coarse = magnus(p0,p1,a0,final,L,mode,64)
        errors["magnus_refinement"] = max(errors["magnus_refinement"],difference(e,coarse))
        d = multiply(e,add(identity(),d1,a0))
        errors["evolution_residual"] = max(errors["evolution_residual"],
            difference(row["evolution"],e),difference(row["fragmentation"],d))
        c1 = add(d1,p0,L)
        c2 = add(scaled(p1,L),multiply(p0,d1),L)
        c2 = add(c2,multiply(p0,p0),L*L/2)
        if mode != "fixed":
            c2 = add(c2,p0,L*L/3)
        errors["coefficient_residual"] = max(errors["coefficient_residual"],
            difference(row["coefficient1"],c1),difference(row["coefficient2"],c2))
        if row["step_refinement_shift"] > result["step_tolerance"]:
            raise RuntimeError("C++ step refinement failed")
        if N == 2:
            for key in ("evolution","fragmentation"):
                errors["momentum_residual"] = max(errors["momentum_residual"],
                    max(abs(sum(r)-1) for r in row[key]))
    if errors["coefficient_residual"] > 1e-7 or max(
            value for key,value in errors.items() if key != "coefficient_residual") > 2e-9:
        raise RuntimeError(f"Independent moment reference failed: {errors}")
    result["independent_checks"] = errors
    result["provenance"] = {
        "source_sha256":source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256":hashlib.sha256(executable.read_bytes()).hexdigest(),
        "python_version":platform.python_version(),"platform":platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(errors,indent=2))
    print("Independent coupled NLO moment and running checks passed")


if __name__ == "__main__":
    main()
