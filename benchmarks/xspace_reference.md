# Coupled x-space reference

Version 0.9.0, AppleClang 21, 53-bit Real. The input is `configs/reference.conf`;
its alpha is interpreted here as the **MS̄ input at mu_0**, without a coupling
conversion. The model has one unit-charge lepton family and two-loop running.
The [contract](../docs/xspace_evolution.md) describes energy weighting,
endpoint projection, comparison bounds and the FF-only observable scope.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel 4
./build/leptonff_xspace configs/reference.conf
python3 python/check_xspace_reference.py build/leptonff_xspace configs/reference.conf \
  --output results/xspace_reference.json
```

The optional driver argument sets the baseline bin count; it must be a power
of two <=32768 so that all refinements fit. The independent CTest uses 64.
Every run also emits a 64-cell snapshot containing all weighted kernels,
boundaries, propagators, FFs and first/second expansion coefficients.
The Python manifest records source/executable/config hashes and environment.

| Grid | Span | FF pair rate after O(alpha) FF matching | Max N=5 moment residual | Max N=8 residual |
| ---: | ---: | ---: | ---: | ---: |
| 4096 | 40 | 0.796732596347 | 1.245e-5 | 5.703e-5 |
| 8192 | 40 | 0.796737290056 | 3.634e-6 | 1.684e-5 |
| 16384 | 40 | 0.796738959018 | 1.045e-6 | 4.884e-6 |
| 32768 | 80 | 0.796738959018 | 1.045e-6 | 4.884e-6 |

The last grid refinement shifts the matched FF rate by `1.669e-6`. Doubling
the span at fixed spacing changes it by about `7e-16`. These diagnose different
effects. The N=2 and N=3 moment residuals against the continuum matrix are below
`2e-13`; N=2 is fixed by the projection identity and N=3 tests its linear-z
moment preservation. Higher moments supply a nontrivial convergence check.

At 16384 cells/span 40 the weighted wrap comparison is `1.10e-13`, the
weighted tail comparison is `9.53e-19`, and imaginary weight mass is `2.74e-15`.
At span 80 the respective wrap/tail comparisons fall to `1.03e-26` and
`7.48e-36`. Signed negative weight mass is about `4.65e-10` and is retained.
These quantities are not the total numerical error or an uncertainty on a
physical fiducial NLL prediction.

| Independent 64-cell comparison | Maximum absolute residual/shift |
| --- | ---: |
| Cell quadrature vs C++ | 8.37e-11 |
| Independent cell refinement 24→40 Gauss nodes | 4.40e-10 |
| Causal power-series propagation vs FFT/RK4 | 9.40e-14 |
| Power-series degree 16→20 | below floating-point resolution |
| First/second lattice expansion coefficients | 3.84e-9 |
| Literal fiducial pair convolution/match | 1.14e-12 |
| Independent continuum O(alpha) FF coefficient | 9.87e-10 |

All **31 CTest checks pass**, including the earlier moment, LL and hard-decay
regressions. Additional checks cover signed endpoint action, invalid domains,
zero coupling, higher-moment refinement and quadratic onset of the FF-matched
remainder. This benchmark includes no H1: the next stage must reconstruct
hard photon-parent terms and align coupling/pair selection before claiming
full fiducial NLL accuracy.
