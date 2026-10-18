# Single-muon spectrum reference

This is the archived v0.13 benchmark and its original manifest. Version 0.14
addresses the panel sensitivity below; see the
[current benchmark](endpoint_quadrature_reference.md). These reproduction
commands describe the v0.13 source, not the current executable.

Version 0.13.0 · Release / AppleClang 21 · 53-bit Real · 17 October 2026

This is the one-family, unpolarised pure-vector parent-at-rest reference, with
M=91.1876 GeV and the original pT>25 GeV, |eta|<2.4, 60<m_pair<120 GeV
pair selection. Both muons fill every selected inclusive distinct-parent pair.
The [measurement contract](../docs/lepton_spectra.md) defines the mean one-leg
cumulatives, hard matching and limitations. No NNLO production, physical Z
spin/flavours, exclusive extra-pair selection or absolute cross section is supplied.

## Reproduction

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 2 --scan > results/lepton_spectrum_scan.json
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 8 > results/lepton_spectrum_panels8.json
python3 python/check_lepton_spectrum_reference.py build/leptonff_lepton_spectrum \
  configs/reference.conf configs/lepton_bins.conf \
  --input results/lepton_spectrum_scan.json \
  --panel-input results/lepton_spectrum_panels8.json \
  --output results/lepton_spectrum_reference.json
```

The saved [manifest](../results/lepton_spectrum_reference.json) includes the raw
points, independent coarse fixture, bin shifts and source/executable/configuration
hashes. The [CTest log](../results/lepton_spectrum_ctest.log) records all 46 passing checks.

## Central bin integrals

The scan central point uses 65536 FF bins over log span 80 (the same spacing
as 32768/40), four FF panels, soft/collinear slivers 2.5e-7, and hard absolute/
relative tolerances 5e-7. Integrals are dimensionless; the driver also reports
densities per GeV or per eta unit. Printed digits are reproducibility digits,
not a claim of that numerical accuracy.

| pT bin [GeV] | Matched | Full O(alpha) |
| --- | --- | --- |
| 25 to 30 | 0.166035461317 | 0.165894273276 |
| 30 to 35 | 0.227208010467 | 0.227527105740 |
| 35 to 40 | 0.323348838884 | 0.323983168055 |
| 40 to 45 | 0.619361560228 | 0.620390931346 |
| 45 to 50 | 0.277847572884 | 0.275398744663 |
| **Histogram sum = twice pair rate** | **1.613801443782** | **1.613194223080** |

| Signed eta bin | Matched | Full O(alpha) |
| --- | --- | --- |
| -2.4 to -1.2 | 0.001920565735 | 0.001852501078 |
| -1.2 to -0.8 | 0.157557854892 | 0.157197007962 |
| -0.8 to 0 | 0.647422301264 | 0.647547602499 |
| 0 to 0.8 | 0.647422301264 | 0.647547602499 |
| 0.8 to 1.2 | 0.157557854892 | 0.157197007962 |
| 1.2 to 2.4 | 0.001920565735 | 0.001852501078 |
| **Histogram sum = twice pair rate** | **1.613801443782** | **1.613194223080** |

The corresponding pair rates are 0.806900721891 matched and
0.806597111540 at full O(alpha). The pT endpoint is M/2=45.5938 GeV;
the 45–50 GeV bin contains the endpoint peak. Signed eta reflection is exact
in this model. All additive hard/FF components and C1/C2 coefficients obey
the same two-entry normalization at every scale point. Couplings and numerical
diagnostics are shared metadata.

## Bin-specific numerical changes

| Refinement | Max absolute pT bin shift | Max absolute eta bin shift | Pair-rate shift |
| --- | --- | --- | --- |
| FF panels 2 to 4 at 8192 bins | 5.683e-05 | 1.431e-06 | +2.333e-08 |
| 8192 to 16384, four panels | 5.482e-06 | 1.531e-06 | +5.098e-07 |
| 16384 to 32768, four panels | 3.329e-05 | 1.990e-07 | +2.935e-07 |
| Span 40 to 80 at fixed spacing | 1.776e-15 | 1.221e-15 | +9.992e-16 |
| Hard slivers/tolerances divided by four | 6.777e-09 | 3.011e-09 | -2.196e-09 |
| FF panels 4 to 8 at 8192 bins | 1.164e-04 | 5.791e-07 | +5.564e-08 |

The additional eight-panel point is a panel audit on the baseline grid, not
a replacement for the central point with its finer grid and hard controls.
Its 4-to-8 panel change is also nonmonotonic: the largest pT shift is 1.16e-4
absolute, and the 45–50 GeV shift is 1.01e-4 (about 3.65e-4 relative).
Stable endpoint FF quadrature remains the leading numerical follow-up; this
study does not certify a smaller pT error from the much smaller total shift.
The pT grid change is nonmonotonic and concentrated across the 45 GeV edge.
Its last shift is about 1.20e-4 relative to the 45–50 GeV integral, while the
pair-rate shift is only 2.94e-7 absolute. Eta has a smaller absolute shift,
but the outer bin is small: 1.99e-7 is about 1.04e-4 relative there. The
conservative mass cells do not smear pT thresholds; angular pT contacts remain
grid sensitive at higher orders after restoring the continuum O(alpha) term.
These are measured shifts, not rigorous error bounds or an extrapolated
continuum value. Finer or narrower endpoint bins need their own convergence study.

## Correlated scale changes

Each row is one variation of the complete spectrum. Entries are absolute
bin-integral shifts from the scan central point. Eta negative bins have
identical shifts to their reflected positive bins.

| Variation | pT 25–30 | 30–35 | 35–40 | 40–45 | 45–50 | Pair rate |
| --- | --- | --- | --- | --- | --- | --- |
| scale_half | +1.864e-05 | -1.887e-05 | -3.543e-05 | -7.501e-05 | +1.402e-04 | +1.474e-05 |
| scale_double | -2.052e-05 | +2.791e-05 | +4.960e-05 | +8.226e-05 | -2.011e-04 | -3.094e-05 |
| start_half | +5.397e-06 | -9.311e-06 | -2.676e-05 | -1.148e-04 | +1.515e-04 | +3.018e-06 |
| start_double | -4.820e-06 | +1.895e-05 | +4.293e-05 | +1.378e-04 | -2.085e-04 | -6.810e-06 |

| Variation | eta 0–0.8 | 0.8–1.2 | 1.2–2.4 |
| --- | --- | --- | --- |
| scale_half | -1.477e-05 | +2.095e-05 | +8.560e-06 |
| scale_double | +1.338e-05 | -3.416e-05 | -1.017e-05 |
| start_half | -2.216e-06 | +1.092e-06 | +4.142e-06 |
| start_double | +4.987e-06 | -6.402e-06 | -5.395e-06 |

Fragmentation and starting scales vary by 1/2 and 2; the hard coupling stays
at M. For starting-scale variations the MSbar input is transported first,
and independent analytic running checks that transport. Opposite bin shifts
cancel substantially in the total. These bands omit unknown constant
H2/d2 contributions, physical flavours/thresholds, spin/current, boosts,
off-shell production and host uncertainties; they are not a full physical
uncertainty budget or evidence for per-mille phenomenology.

## Verification

C++ tests check the one-lepton mean (including pairs where only one muon
clears the extra threshold), analytic Born/Jacobian bins, a signed literal
non-endpoint tensor, asymmetric eta bins crossing zero, endpoint contacts,
reflection and two-entry closure. Per-bin first-order scale cancellation and
absence of a retained alpha² mass-log scale shift pass. Representative pT
and eta cumulatives agree with the massive real/virtual reference within
5e-6 in rate, using separately controlled hard integration.

Python independently builds causal FFs and z-quadrature hard source cells,
evaluates literal signed lattice pairs and mass-cell intersections, and uses
direct allowed-azimuth integration for selected muon/photon geometry samples.
The leading-power massive spectrum plus that geometry gives a separate
continuum O(alpha) bin reference. It samples the full resummed geometry
matrices; inherited coupled contraction tests and bin refinements cover the
remaining action. No histogram renormalization or signed-weight clipping is used.

| Independent coarse check | Maximum absolute residual |
| --- | --- |
| Causal FF weights | 1.514e-13 |
| Hard source cells | 5.077e-09 |
| Literal FF bin action | 1.023e-13 |
| Muon hard projection | 1.839e-09 |
| Photon hard projection | 6.200e-13 |
| Selected hard-geometry coefficient samples | 2.019e-05 |
| Allowed-azimuth sample refinement | 6.072e-05 |
| Continuum O(alpha) bin coefficient | 2.099e-06 |
| Continuum angular refinement | 3.071e-06 |
| Full O(alpha) bin integral | 2.437e-09 |

The signed literal tensor guards non-endpoint lattice indexing, which a
Born-only or telescoping closure test cannot establish. Explicit angular
contacts cover both energy branches, including narrow regions where only the
higher-energy muon passes the threshold. Fixed-order and dressed-hard endpoint
geometry agree. Nearly coincident generated contacts are merged below
quadrature resolution; physical integration endpoints and source-distribution
actions retain their original definitions.

Physical flavours/spin/boosts, a full uncertainty budget, an external aligned
comparison and a validated NNLO production adapter remain outstanding.
