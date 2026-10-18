# Matched bare-pair mass spectrum

Version 0.12.0, Release, AppleClang 21, 53-bit Real; **43 CTest checks pass**.
One-family MSbar model, unpolarised vector at rest, inclusive opposite-charge
pairs from distinct hard directions. Normalisation is massless inclusive
vector Born. This is not a physical off-shell DY spectrum.

M=91.1876 GeV, m=mu_0=0.1056583755 GeV,
alpha_MSbar(mu_0)=0.0072973525692838015, pT>25 GeV, |eta|<2.4,
60<m_pair<120 GeV. No on-shell coupling conversion, physical thresholds or
extra-pair veto. See the [measurement contract](../docs/mass_spectrum.md).

```sh
./build/leptonff_mass_spectrum configs/reference.conf configs/mass_bins.conf \
  8192 4 --scan > results/mass_spectrum_fine_scan.json
python3 python/check_mass_spectrum_reference.py build/leptonff_mass_spectrum \
  configs/reference.conf configs/mass_bins.conf \
  --input results/mass_spectrum_fine_scan.json --output results/mass_spectrum_reference.json
```

Checked JSON includes source/executable/configuration/bin-config/input hashes,
the independent fixture, components and correlated numerical/scale shifts.
Integrals are dimensionless; density fields divide by the width in GeV.

## Central spectrum

Final point: 65536 FF bins over log span 80, eight FF panels, soft/collinear
slivers 2.5e-7, geometry tolerances 5e-7, mu_frag=M and hard alpha_MSbar(M).
Resolved total-loss cells are integrated conservatively; zero loss stays an atom.

| Pair mass [GeV] | Matched bin integral | Full O(alpha) bin integral |
| --- | --- | --- |
| 60–70 | 0.001901525769 | 0.001697161225 |
| 70–80 | 0.014582084994 | 0.014718868692 |
| 80–85 | 0.019338178060 | 0.020312960814 |
| 85–90 | 0.064351528364 | 0.071346121806 |
| 90–92 | 0.706727440407 | 0.698521999004 |
| 92–100 | 0 | 0 |
| 100–120 | 0 | 0 |
| **Sum / fiducial rate** | **0.806900757593** | **0.806597111540** |

The pair mass cannot exceed fixed M. Empty upper bins are a model property,
not a hadronic Z line shape. An interior edge at M assigns Born and its soft
neighbourhood to the lower bin. Digits describe reproducibility, not accuracy.

## Bin-specific numerical convergence

| FF bins / span | Panels | Total | 85–90 GeV | 90–92 GeV |
| --- | --- | --- | --- | --- |
| 4096 / 40 | 4 | 0.806899290812 | 0.064256660745 | 0.706820260375 |
| 4096 / 40 | 8 | 0.806899382807 | 0.064256532103 | 0.706820342249 |
| 8192 / 40 | 8 | 0.806899976382 | 0.064329586322 | 0.706748811330 |
| 16384 / 40 | 8 | 0.806900477402 | 0.064344843425 | 0.706733800796 |
| 32768 / 40 | 8 | 0.806900756447 | 0.064351528239 | 0.706727440342 |
| 65536 / 80 | 8 | 0.806900756447 | 0.064351528239 | 0.706727440342 |
| 65536 / 80, refined hard controls | 8 | 0.806900757593 | 0.064351528364 | 0.706727440407 |

Point assignment initially changed the 85–90 GeV bin by 2.79e-4 in an
8192-to-16384 comparison despite a total shift near 2.5e-7. Cell integration
reduces that bin shift to 1.53e-5, then 6.69e-6 in the next doubling.

| Last refinement in final scan | Maximum absolute bin shift | Total shift |
| --- | --- | --- |
| FF panels 4 -> 8 at 8192 bins | 1.61e-7 | +5.56e-8 |
| Grid 8192 -> 16384 | 1.53e-5 | +5.01e-7 |
| Grid 16384 -> 32768 | 6.69e-6 | +2.79e-7 |
| Span 40 -> 80, fixed spacing | 6.67e-16 | +1.00e-15 |
| Hard slivers/tolerances divided by four | 9.32e-10 | +1.15e-9 |

The last grid shift is about 1.04e-4 relative to the 85–90 GeV bin and dominates
its measured numerical changes. These are convergence diagnostics, not rigorous
bounds or evidence for twelve digits. Narrower bins need their own study.
Finite-grid totals can differ from the older nodal rate driver; both converge
to the same observable and restore the same continuum O(alpha) target.

## Correlated scale shifts

Absolute bin-integral shifts from the final central point; each row is a
single common variation of the full spectrum:

| Variation | 60–70 | 70–80 | 80–85 | 85–90 | 90–92 | Total |
| --- | --- | --- | --- | --- | --- | --- |
| mu_frag/M=0.5 | +1.679e-5 | +1.920e-5 | -2.592e-5 | -5.145e-4 | +5.192e-4 | +1.473e-5 |
| mu_frag/M=2 | -2.201e-5 | -9.212e-6 | +5.896e-5 | +7.004e-4 | -7.590e-4 | -3.094e-5 |
| mu_0 / nominal=0.5 | +6.139e-7 | +8.813e-6 | -3.427e-7 | -1.858e-4 | +1.797e-4 | +3.018e-6 |
| mu_0 / nominal=2 | -1.908e-7 | -2.628e-6 | +1.638e-5 | +2.886e-4 | -3.089e-4 | -6.810e-6 |

Hard alpha stays at M. For mu_0 changes the MSbar input is transported from
the nominal scale before recomputing boundaries/evolution; independent analytic
running verifies this. No physical threshold prescription is implied.
Large opposite endpoint-bin shifts cancel in the total. Its narrow band is
therefore not each bin's uncertainty. These are model diagnostics, not a full
uncertainty estimate: unknown H2/d2 constants, physical flavours/spin/current,
off-shell production and host errors remain absent. No physical per-mille claim.

## Verification

All additive components and first/second coefficients sum to varied totals.
Independent coarse component/expansion closure is within 2.85e-13; the final
matched closure is zero at printed precision. Couplings/diagnostics are shared
metadata. No signed contribution is clipped and no bins are renormalised.

| Independent coarse check | Maximum absolute residual |
| --- | --- |
| Causal-series FF weights | 1.52e-13 |
| Hard source cells | 5.08e-9 |
| Literal bin-integrated FF action | 2.32e-13 |
| Muon hard mass projection | 3.55e-9 |
| Photon hard mass projection, both charges | 2.15e-13 |
| Continuum O(alpha) bin coefficient | 1.43e-6 |
| Continuum angular refinement | 7.47e-7 |
| Full O(alpha) bin integral | 1.66e-9 |

The independent script uses z-cell quadrature, causal series and literal cell
intersections; production uses FFT/RK evolution and cumulative/prefix sums.
The O(alpha) reference uses the leading-power massive spectrum plus direct
allowed-azimuth geometry. It does not independently compute every resummed
geometry bin. Inherited current/fixed-fraction/total-rate checks, direct-tail
tests, bin-specific refinements and per-bin mass-log scale checks supply the
remaining evidence within this restricted theory.

The subsequent [v0.13 benchmark](lepton_spectrum_reference.md) adds single-muon
pT/eta output. A full uncertainty budget, physical flavours, spin/boosts and
an NNLO production adapter remain outstanding. This v0.12 manifest records
its original source revision; current regression tests also exercise it.
