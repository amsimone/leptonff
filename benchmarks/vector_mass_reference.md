# Angular-inclusive vector hard-matching reference

Reproduce with:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
python3 python/check_vector_mass_reference.py build/leptonff_hard configs/vector_mass.conf \
  --output results/vector_mass_reference.json
```

Version 0.4.0, AppleClang 21, 53-bit `long double`; inputs are in
[vector_mass.conf](../configs/vector_mass.conf). The JSON records compiler,
integration controls, all generated configurations, scale/mass scans, independent
checks, and SHA-256 hashes of sources, executable and base configuration.

All rates below are relative to massless vector Born, with fixed alpha and
M=91.1876 GeV. The first three rows use m_mu=0.1056583755 GeV.

| Mass window / fixture | Leading-power hard+FF O(alpha) | Exact massive O(alpha) |
| --- | ---: | ---: |
| 60–120 GeV | 0.986414026151 | 0.986414112746 |
| 70–85 GeV, real radiation only | 0.051195075920 | 0.051195019127 |
| 0–120 GeV, inclusive | 1.001742114599 | 1.001742142656 |
| 60–120 GeV, m_mu=2 GeV | 0.994719175084 | 0.994737499211 |

At the central scale mu_frag=M, the 60–120 GeV coefficients multiplying
alpha/(2*pi) are H1=4.236705026008 and 2*D1=-15.934535227315. Their sum is
-11.697830201307; the inclusive coefficient is exactly 3/2. Varying mu_frag/M
through {0.5,1,2} changes the central mass-window rate by 6.09e-13, while the
separate coefficients vary substantially. This tests an O(alpha) cancellation,
not an uncertainty estimate for omitted higher orders.

| Muon mass [GeV] | Massive minus leading-power rate, 60–120 GeV |
| ---: | ---: |
| 0.1056583755 | 8.660e-8 |
| 0.05282918775 | 2.352e-8 |
| 0.026414593875 | 6.350e-9 |
| 0.0132072969375 | 1.705e-9 |

Tightening the central C++ tolerances by ten shifts the leading-power rate by
3.55e-13 and the reported massive rate by less than the displayed precision.
The independent Python coefficient residual is below 3.5e-10 for leading power
and 8.9e-8 for the inclusive/massive checks, after extrapolating soft separators
1e-4 and 5e-5. The real-tail coefficient residual is below 5e-13. These checks
use distinct quadratures and virtual formulae; they are numerical evidence for
this fixture, not guarantees for arbitrary endpoint bins.

The C++ suite additionally integrates a separate Dirac-trace matrix element
and compares it with the exact scalar radiator. The source discrepancy found
during this check and the equations retained are recorded in
[vector matching](../docs/vector_mass_matching.md).

This benchmark has no angular cuts, axial current, resummation or hadronic
production. It does not close the full fiducial matched-decay gate G2.
