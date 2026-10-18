# NLO time-like charged-lepton kernel reference

Version 0.7.0, AppleClang 21, 53-bit Real. This benchmark checks the one-family,
unit-charge MS̄ lepton block. Values are coefficients of `(alpha/(2*pi))²`,
not probabilities or evolved fragmentation functions. The
[source audit](../docs/nlo_kernel_audit.md) pins formulas and charge mapping.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel 4
./build/leptonff_kernels
python3 python/check_nlo_kernel_reference.py build/leptonff_kernels \
  --output results/nlo_kernel_reference.json
```

The independent script records source/executable hashes and environment
metadata. Its checks use integral-defined dilogarithms, transformed fixed
Gauss quadrature and analytic harmonic-sum moments. Production C++ uses a
dilogarithm series/reflections and adaptive Gauss–Kronrod integration.

| Kernel/component | N=1 | N=2 | N=3 |
| --- | ---: | ---: | ---: |
| Photonic G | −0.719311605774 | 2.794410824823 | 3.262199136479 |
| Non-singlet pair B | 0 | 1.185185185185 | 1.921296296296 |
| Singlet pair S | divergent | −0.592592592593 | −0.203287037037 |
| Same-charge interference I | 1.438623211547 | 0.681761996250 | 0.488666451533 |
| Opposite-charge interference J | 0.719311605774 | 0.037549609525 | 0.005465956870 |
| Same charge G+B+S+I | divergent | 4.068765413665 | 5.468874847270 |
| Opposite charge S+J | divergent | −0.555042983068 | −0.197821080167 |
| Valence G+B+I−J | 0 | 4.623808396733 | 5.666695927437 |

The raw valence N=1 residual is `7.07e-14`. The source-separated photonic
first moment is nonzero; the conserved charge check uses the complete valence
combination. The singlet is cancelled analytically before that action.
The C++ API rejects divergent N=1 inclusive-lepton moments rather than
subtracting divergent numbers. The output represents them as JSON `null`.

The diagnostic also checks N=5 and each kernel's finite-bound action
`integral_[0.37,1] dz z² P1(z)`, including delta and plus terms. Selected values
are `5.615681759580` for the same-charge channel and `5.766426409939` for
valence. This exposes the lower-bound endpoint subtraction separately from
full-interval moments.

| Numerical check | Maximum absolute difference |
| --- | ---: |
| C++ moment tolerance refinement 1e-10 → 1e-11 | 3.15e-11 |
| Independent Gauss refinement 48 → 80 nodes | 2.86e-11 |
| Independent integral vs C++ | 2.10e-12 |
| Analytic moment checks | 2.10e-12 |

All 24 CTest checks pass, including the earlier decay/LL/matching references.
The new block is not connected to those LL solvers. NLO photon channels,
per-parent momentum conservation, running evolution and observable pair
selection remain required before a coupled NLL prediction.

This records the version 0.7 lepton-block milestone. Version 0.8 extends the
same driver/checker commands with the photon entries and supplies coupled
running real moments; see the [current benchmark](moment_reference.md).
