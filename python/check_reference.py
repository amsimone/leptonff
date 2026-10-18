#!/usr/bin/env python3
"""Independent Simpson/endpoint-transformed check of the reference geometry.

No third-party modules. This is not a general fiducial integrator: it requires a
parent inside the mass window and an eta cut looser than the Born pt bound.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess


def simpson(function, lower, upper, tolerance=2e-10):
    def recurse(lo, hi, f_lo, f_mid, f_hi, whole, tol, depth):
        mid = (lo + hi) / 2
        left_mid = function((lo + mid) / 2)
        right_mid = function((mid + hi) / 2)
        left = (mid - lo) * (f_lo + 4 * left_mid + f_mid) / 6
        right = (hi - mid) * (f_mid + 4 * right_mid + f_hi) / 6
        correction = left + right - whole
        if abs(correction) <= 15 * tol:
            return left + right + correction / 15
        if depth == 0:
            raise RuntimeError("Independent Simpson reference did not converge")
        return (
            recurse(lo, mid, f_lo, left_mid, f_mid, left, tol / 2, depth - 1)
            + recurse(mid, hi, f_mid, right_mid, f_hi, right, tol / 2, depth - 1)
        )

    f_lo = function(lower)
    f_mid = function((lower + upper) / 2)
    f_hi = function(upper)
    initial = (upper - lower) * (f_lo + 4 * f_mid + f_hi) / 6
    return recurse(lower, upper, f_lo, f_mid, f_hi, initial, tolerance, 40)


def independent_coefficient(result):
    mass = result["parent_mass"]
    ratio = 2 * result["pt_min"] / mass
    born = math.sqrt(1 - ratio * ratio)
    if not (0 < ratio < 1 and result["mass_min"] < mass < result["mass_max"]
            and math.tanh(result["eta_max"]) > born):
        raise ValueError("Configuration outside the independent reference geometry")
    lower = max(ratio, (result["mass_min"] / mass) ** 2)
    log_scale = 2 * math.log(result["mu_frag"] / result["muon_mass"])

    # Below the cut: F(z)=0. Above: rationalise A(z)-A(1) before multiplying
    # the singular kernel. The sin^2 map resolves the square-root threshold
    # and soft endpoint, unlike the production Gauss-Kronrod implementation.
    below = simpson(
        lambda z: -born * (1 + z*z) / (1-z)
        * (log_scale - 1 - 2*math.log1p(-z)), 0, lower
    )

    def transformed(u):
        if u == 0 or u == 1:
            return 0.0
        angle = math.pi * u / 2
        sine, cosine = math.sin(angle), math.cos(angle)
        one_minus_z = (1-lower) * cosine*cosine
        z = 1 - one_minus_z
        acceptance = math.sqrt(max(0.0, 1 - (ratio/z)**2))
        jacobian = (1-lower) * math.pi * sine * cosine
        rationalised = -(1+z*z) * ratio*ratio * (1+z) / (z*z*(acceptance+born))
        return rationalised * (log_scale - 1 - 2*math.log(one_minus_z)) * jacobian

    return born, 2 * (below + simpson(transformed, 0, 1))


def source_digest(root):
    paths = [root / "CMakeLists.txt"]
    for directory in ("include", "src", "apps", "python", "tests"):
        paths.extend(path for path in (root / directory).rglob("*")
                     if path.is_file() and path.suffix in (".hpp", ".cpp", ".py"))
    digest = hashlib.sha256()
    for path in sorted(paths):
        digest.update(str(path.relative_to(root)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("config", type=Path)
    parser.add_argument("--output", type=Path, help="Write results plus provenance JSON")
    args = parser.parse_args()
    executable, config = args.executable.resolve(), args.config.resolve()
    run = subprocess.run([str(executable), str(config)], check=True,
                         capture_output=True, text=True)
    result = json.loads(run.stdout)
    born, coefficient = independent_coefficient(result)
    residual = coefficient - result["ff_coefficient_of_alpha_over_2pi"]
    if abs(born - result["born_acceptance"]) > 1e-12 or abs(residual) > 1e-7:
        raise RuntimeError(f"Independent reference mismatch: delta coefficient={residual}")
    if abs(result["lepton_momentum_coefficient_residual"]) > 1e-6:
        raise RuntimeError("Lepton momentum sum rule failed")
    if abs(result["photon_momentum_coefficient_residual"]) > 1e-6:
        raise RuntimeError("Photon momentum sum rule failed")
    result["independent_coefficient"] = coefficient
    result["independent_coefficient_residual"] = residual
    result["provenance"] = {
        "source_sha256": source_digest(Path(__file__).resolve().parents[1]),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"Born acceptance             {born:.12f}")
    print(f"O(alpha/2pi) FF coefficient  {coefficient:.12f}")
    print(f"Born plus O(alpha) FF       {result['born_plus_ff_order1']:.12f}")
    print(f"Independent coefficient Δ  {residual:.3e}")
    print("Independent reference and momentum checks passed")


if __name__ == "__main__":
    main()
