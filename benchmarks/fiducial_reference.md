# Fiducial vector hard+FF and additive LL benchmark

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
python3 python/check_fiducial_reference.py build/leptonff_fiducial configs/reference.conf \
  --output results/fiducial_reference.json
```

Version 0.5.0, AppleClang 21, 53-bit Real. The input fixture is the existing
[reference.conf](../configs/reference.conf): M=91.1876 GeV,
m_mu=0.1056583755 GeV, pT>25 GeV, |eta|<2.4, 60<m_mumu<120 GeV.
All rates are dimensionless relative to total massless vector Born.

| Quantity | Rate |
| --- | ---: |
| Fiducial Born | 0.836268494675 |
| FF-only O(alpha), previous driver | 0.796111391223 |
| Full leading-power O(alpha), hard+FF | 0.806597111174 |
| Exact massive O(alpha), finite soft sliver omitted | 0.806595873127 |
| Additive NLO+LL, baseline lattice | 0.807231636355 |
| Additive NLO+LL, refined grid | 0.807231627694 |
| Additive NLO+LL, refined soft cutoff | 0.807232541918 |

At mu_frag=M, the coefficients of alpha/(2*pi) are:

| Piece | Coefficient |
| --- | ---: |
| Projected mass-spectrum hard coefficient | 7.796947697227 |
| Finite geometry remainder | 1.231493906744 |
| Total fiducial hard coefficient | 9.028441603971 |
| Two FF insertions | -34.576172659480 |
| Full leading-power coefficient | -25.547731055509 |

Varying mu_frag/M through {0.5,1,2} changes the **fixed-order** rate by
6.51e-13. The higher-order resummed result has no corresponding cancellation
claim or NLL scale interpretation.

The finite remainder starts with epsilon=delta=1e-6 and quadrature absolute/
relative tolerances 2e-6. The final result halves both slivers and quarters the
tolerances. Separate rate shifts are 1.79e-10 (soft), 4.49e-10 (collinear) and
6.25e-10 (quadrature). The massive reference shift is -9.64e-9. The LL grid
shift is -8.66e-9 and its soft-cutoff shift is +9.14e-7; these are distinct
numerical controls. Default lattice settings are unchanged from the FF benchmark.

| Muon mass [GeV] | Massive minus leading-power rate |
| ---: | ---: |
| 0.1056583755 | -1.23805e-6 |
| 0.05282918775 | -3.08546e-7 |
| 0.026414593875 | -6.93079e-8 |

The mass scan approaches leading power. Its smallest differences also depend
on the finite soft sliver; the manifest reports the controls explicitly.

Python uses independent Gauss-Legendre rules and direct azimuthal integration,
rather than the production spherical-cap formula and adaptive GK quadrature.
Its 48/96-point geometry estimates differ by 9.73e-7 in the coefficient;
the fine estimate differs from C++ by 1.77e-7. The independent full fixed-order
coefficient agrees within the same amount. This does not establish identical
precision for arbitrary cuts.

C++ checks also cover massless/massive Dirac traces, eta-only cuts, real-tail
bins, the angular-inclusive reduction, mass and regulator scans, the alpha
expansion, no-evolution/zero-coupling limits, and rejection of stale hard data.

The result is an unpolarised pure-vector parent at rest. NLO splitting,
LL dressing of H1, running/mixing, DY spin correlations, boosts and physical
axial mass effects remain outside this benchmark. See the
[derivation and matching contract](../docs/fiducial_vector_matching.md).

A separate eta-only fixture (pT>0, |eta|<0.8, 60–120 GeV) requires a smaller
angular sliver, delta=1e-10. Its massive-minus-leading-power rate differences
are 7.31e-6, 3.64e-6 and 1.81e-6 as the mass is halved twice. Thus sharp tied
angular cuts have linear mass sensitivity; an m_mu^2-only fit is inappropriate.
The broader default-driver slivers are not validated to this accuracy for that
geometry.

Version 0.6 adds the separate [LL-dressed hard benchmark](hard_dressed_reference.md);
the additive numbers above remain a comparison fixture.
