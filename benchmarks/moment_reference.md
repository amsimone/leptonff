# Complete one-family NLO matrix and running-moment reference

Version 0.8.0, AppleClang 21, 53-bit Real. These are MS̄ inclusive fragmentation
moments in the one-family model, not fiducial acceptances. The
[photon derivation](../docs/nlo_photon_kernels.md) and
[evolution contract](../docs/moment_evolution.md) record source conventions,
matrix orientation and expansion content.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel 4
./build/leptonff_kernels
./build/leptonff_moments
python3 python/check_nlo_kernel_reference.py build/leptonff_kernels \
  --output results/nlo_kernel_reference.json
python3 python/check_moment_reference.py build/leptonff_moments \
  --output results/moment_reference.json
```

The saved JSON includes source/executable hashes and the Python environment.
The moment driver reports all nine propagator/FF entries for N=2,3,5 in
fixed, one-loop and two-loop coupling modes, their kernels and boundaries,
both expansion coefficients, running integrals and step-refinement diagnostics.

## Full NLO matrix checks

Photon entries below multiply `(alpha/(2*pi))²`. There are no plus terms in
these entries; photon survival includes its explicit virtual `-delta(1-z)`.

| Channel | N=1 | N=2 | N=3 |
| --- | ---: | ---: | ---: |
| mu→gamma | −7.421533662112 | −3.513722430601 | −2.115378563388 |
| gamma→each mu charge | 1.000000000001 | 1.203703703704 | 1.044120370371 |
| gamma→gamma | divergent | −2.407407407407 | −1.188796296296 |

The existing lepton entries are unchanged. All three parent momentum sums at
N=2 vanish, with maximum residual `3.76e-12`. Independent photon expressions
come from space-like QED plus the full analytic-continuation correction,
including the crossed squared-logarithm term. The complete kernel diagnostic
also checks N=5 and finite-bound actions.

| Kernel check | Maximum absolute difference |
| --- | ---: |
| Independent quadrature vs C++ | 9.69e-12 |
| Independent quadrature 48 → 80 nodes | 2.86e-11 |
| Analytic moment checks | 2.10e-12 |

## Coupled evolution

The input is `alpha_MSbar(mu_0)=1/137.035999084`, `mu_0=m_mu=0.1056583755 GeV`
and `mu_frag=91.1876 GeV`. This defines a model MS̄ coupling input; it is not
an on-shell-to-MS̄ conversion. One lepton family remains active throughout.
The final coupling is `0.00737469262748935` with two-loop running, compared
with `0.00737455598709066` at one loop.

The mu− parent row of the boundary-transported fragmentation matrix is:

| N and running | Observed mu− | Observed mu+ | Observed gamma |
| --- | ---: | ---: | ---: |
| 2, fixed | 0.977257881774 | 0.000043903596 | 0.022698214630 |
| 2, one loop | 0.977151365431 | 0.000044363324 | 0.022804271244 |
| 2, two loops | 0.977151177919 | 0.000044364136 | 0.022804457946 |
| 3, two loops | 0.963624164411 | 0.000012966791 | 0.009082668511 |
| 5, two loops | 0.945917483435 | 0.000003636096 | 0.003935878285 |

The N=2 row sums to one. Number moments of individual inclusive lepton
channels diverge; they are not used as normalisations. Charge-odd evolution
is checked against a separate scalar exponential after analytic singlet
cancellation.

All reference evolutions reached the absolute `1e-11` step-shift criterion
at 32 RK4 steps; the largest actual shift was `3.09e-14`. Kernel quadrature
uses separate `1e-11` absolute/relative tolerances. Independent Python obtains
the coupling from its implicit analytic equation and evolves the matrix in
`log(a/a_initial)` with fourth-order Magnus/exponential steps.

| Independent evolution check | Maximum absolute difference |
| --- | ---: |
| Running coupling and integrated powers | 1.74e-16 |
| Analytic/independent kernel and boundary moments | 3.15e-12 |
| Propagator and boundary-transported FF matrices | 1.13e-14 |
| Magnus refinement 64 → 128 steps | 9.77e-15 |
| First/second expansion coefficients | 1.51e-10 |
| N=2 evolved per-parent momentum | 6.67e-16 |

C++ tests also cover noncommuting matrix composition at a larger coupling,
backward inversion, equal-scale and zero-coupling identities, cubic onset of
the residual after O(α²) subtraction, and invalid/strong-running inputs.
All **28 CTest checks pass**, including the earlier LL and fiducial regressions.

X-space inversion and endpoint control are not supplied by these moments.
Neither the old LL fiducial rates nor their coupling convention are changed.
Fiducial NLL matching still needs hard-coupling alignment, hard-parent channel
counting and a documented extra-pair selection. Physical thresholds, hadronic
running and spin-correlated boosted production remain outside this reference.
