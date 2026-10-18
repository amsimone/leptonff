# Endpoint quadrature reference

Version 0.14.0 · Release / AppleClang 21 · 53-bit Real · 17 October 2026

The [numerical contract](../docs/endpoint_quadrature.md) specifies direct endpoint
geometry and a common FF panel partition aligned with pT/eta contacts. The
physical model and hard/FF accuracy are unchanged: one lepton family, a fixed
unpolarised vector at rest, inclusive distinct-parent opposite-charge pairs,
and no extra-pair veto. Both muons enter every selected pair histogram. This
is not a physical NNLO Drell–Yan prediction.

## Reproduction

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 2 --scan > results/lepton_spectrum_endpoint_scan.json
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 8 > results/lepton_spectrum_endpoint_panels8.json
python3 python/check_lepton_spectrum_reference.py build/leptonff_lepton_spectrum \
  configs/reference.conf configs/lepton_bins.conf \
  --checked-coarse build/lepton_spectrum_endpoint_checked_coarse.json \
  --input results/lepton_spectrum_endpoint_scan.json \
  --panel-input results/lepton_spectrum_endpoint_panels8.json \
  --output results/lepton_spectrum_endpoint_reference.json
```

The [manifest](../results/lepton_spectrum_endpoint_reference.json) records all
production points, the independent coarse fixture and source/executable/input
hashes. CTest writes the independently checked coarse artifact into the build
directory. `--checked-coarse` reuses it only when source, executable and both
configuration hashes match; omit that option to recompute the independent
angular audit. The [test log](../results/endpoint_quadrature_ctest.log) records
all 47 passing tests. The separate [v0.13 manifest](../results/lepton_spectrum_reference.json)
remains an archived calculation at its original source revision.

## Panel stability at fixed FF lattice

All points in this table use 8192 bins over log span 40 and the same hard
controls. The base panel count is augmented by the same physical contacts.

| Base FF panels | Active FF nodes | 40–45 GeV | 45–50 GeV | Pair rate |
| --- | --- | --- | --- | --- |
| 2 | 33 | 0.619315447015 | 0.277883900134 | 0.806900112540 |
| 4 | 41 | 0.619315447151 | 0.277883900141 | 0.806900099768 |
| 8 | 53 | 0.619315445486 | 0.277883900692 | 0.806900088912 |

The maximum four-to-eight-panel pT bin change is **3.245e-08**,
compared with **1.164e-04** in v0.13: a reduction by about 3586.
No signed contribution is clipped and no bins are renormalised. The reduction
comes from preserving forbidden pT support in each interpolation panel and
evaluating the endpoint strips directly.

The largest C2 change under these panel refinements is 0.000e+00.
Its hard geometry now uses the original lattice endpoint strips. The remaining
compressed resolved/resolved geometry first enters the matched rate at
O(alpha³). C2 can still change under lattice or hard-integration refinement.

## Grid, domain and hard-integration study

| Refinement | Max absolute pT bin shift | Max absolute eta bin shift | Pair-rate shift |
| --- | --- | --- | --- |
| Base FF panels 2 to 4 | 2.940e-08 | 1.278e-08 | -1.277e-08 |
| 8192 to 16384, four panels | 2.465e-05 | 1.635e-06 | +4.492e-07 |
| 16384 to 32768, four panels | 1.114e-05 | 2.248e-07 | +2.578e-07 |
| Span 40 to 80, fixed spacing | 1.554e-15 | 1.221e-15 | +9.992e-16 |
| Hard slivers/tolerances divided by four | 3.491e-09 | 8.830e-10 | +1.316e-09 |
| Base FF panels 4 to 8 | 3.245e-08 | 5.822e-09 | -1.086e-08 |

The extra eight-panel audit is on the 8192-bin lattice. It does not replace
a joint finest-grid/panel study for a tighter numerical budget.
The panel interpolation problem is controlled at the tested resolution; the pT
lattice-spacing change is larger. These are measured shifts, not rigorous
error bounds or an extrapolated continuum value. Angular pT thresholds remain
unsmeared, and restoring continuum O(alpha) does not remove higher-order
finite-grid effects. Narrower bins and a tighter endpoint precision target
require their own lattice study. A small total-rate change does not bound
individual bins.

## Central bin integrals

The scan central point uses 65536 bins over log span 80, four base FF panels
plus the common contacts, and soft/collinear cutoffs 2.5e-7 with hard absolute/
relative tolerances 5e-7. It has the same lattice spacing as 32768/40. Integrals
are dimensionless relative to the massless inclusive vector Born rate; density
fields divide by bin width. Printed digits describe reproducibility, not accuracy.

| pT bin [GeV] | Matched | Full O(alpha) |
| --- | --- | --- |
| 25 to 30 | 0.166032594986 | 0.165894273276 |
| 30 to 35 | 0.227212351174 | 0.227527105740 |
| 35 to 40 | 0.323359754660 | 0.323983168055 |
| 40 to 45 | 0.619326942455 | 0.620390931346 |
| 45 to 50 | 0.277869973029 | 0.275398744663 |
| **Sum = twice pair rate** | **1.613801616303** | **1.613194223080** |

| Signed eta bin | Matched | Full O(alpha) |
| --- | --- | --- |
| -2.4 to -1.2 | 0.001920545248 | 0.001852501078 |
| -1.2 to -0.8 | 0.157557419939 | 0.157197007962 |
| -0.8 to 0 | 0.647422842965 | 0.647547602499 |
| 0 to 0.8 | 0.647422842965 | 0.647547602499 |
| 0.8 to 1.2 | 0.157557419939 | 0.157197007962 |
| 1.2 to 2.4 | 0.001920545248 | 0.001852501078 |
| **Sum = twice pair rate** | **1.613801616303** | **1.613194223080** |

The central pair rate is 0.806900808151; the fixed-order target is
0.806597111540. The parent mass is fixed at 91.1876 GeV, so
the pT endpoint is 45.5938 GeV. Signed eta reflection is exact only in this
unpolarised rest-frame model.

## Correlated scale variations

Each row is one variation of the entire spectrum from the scan central point.
Eta negative bins have the same shifts as their reflected positive bins.

| Variation | pT 25–30 | 30–35 | 35–40 | 40–45 | 45–50 | Pair rate |
| --- | --- | --- | --- | --- | --- | --- |
| scale_half | +1.890e-05 | -1.927e-05 | -3.626e-05 | -7.251e-05 | +1.386e-04 | +1.473e-05 |
| scale_double | -2.077e-05 | +2.830e-05 | +5.042e-05 | +7.980e-05 | -1.996e-04 | -3.094e-05 |
| start_half | +5.396e-06 | -9.309e-06 | -2.677e-05 | -1.149e-04 | +1.516e-04 | +3.019e-06 |
| start_double | -4.820e-06 | +1.895e-05 | +4.295e-05 | +1.378e-04 | -2.085e-04 | -6.810e-06 |

| Variation | eta 0–0.8 | 0.8–1.2 | 1.2–2.4 |
| --- | --- | --- | --- |
| scale_half | -1.482e-05 | +2.099e-05 | +8.562e-06 |
| scale_double | +1.343e-05 | -3.420e-05 | -1.017e-05 |
| start_half | -2.216e-06 | +1.093e-06 | +4.142e-06 |
| start_double | +4.987e-06 | -6.402e-06 | -5.396e-06 |

Fragmentation and starting scales vary by 1/2 and 2. The hard coupling remains
at M; starting-scale variations transport the common MSbar input first. Every
additive component and C1/C2 still sums to twice the varied pair total. The
bands omit unknown constant H2/d2 terms, physical flavours and thresholds,
spin/current, boosts, production PDFs and host errors. They are not a complete
physical uncertainty budget or a per-mille accuracy claim.

## Independent verification

The new C++ test uses a lattice where every resolved fraction is too small
for a 45 GeV muon. It compares the full signed muon/photon tensor with literal
lattice evaluation and checks zero forbidden support, exact contact ownership,
panel-independent C2 and the cubic coupling remainder. Existing single-muon
closure, mass-log, scale and massive-reference tests pass.

The independent Python audit reconstructs causal source arrays, nonuniform
cubic projection and mixed geometry contractions. It checks every coarse
endpoint strip with direct allowed-azimuth integration, including both photon
orientations, and reconstructs C1/C2. Its angular integration supplies both
collinear branch changes and splits spherical-cap contact points. Interior
matrices are sampled independently; the remaining interior contraction uses
the separately checked matrix data. This does not independently angularly
reconstruct every production matrix entry.

| Independent coarse check | Maximum absolute residual |
| --- | --- |
| Causal FF weights | 1.514e-13 |
| Hard source cells | 5.077e-09 |
| Literal FF bin action | 1.023e-13 |
| Muon hard projection | 1.839e-09 |
| Photon hard projection | 6.200e-13 |
| Direct endpoint geometry coefficients | 1.258e-05 |
| Endpoint angular refinement | 7.924e-05 |
| Interior/endpoint geometry samples | 2.286e-05 |
| Sample angular refinement | 1.858e-04 |
| Forbidden interpolation support | 0.000e+00 |
| Signed cubic moments | 2.220e-16 |
| Mixed muon geometry algebra | 6.608e-13 |
| Mixed photon geometry algebra | 5.076e-16 |
| Lattice C1 | 1.949e-09 |
| C2 source/endpoint algebra | 1.757e-07 |
| C2 from independent endpoint angles | 2.904e-03 |
| C2 angular refinement | 2.398e-02 |
| Continuum O(alpha) coefficient | 1.383e-06 |
| Continuum angular refinement | 2.547e-06 |
| Full O(alpha) bin rate | 1.607e-09 |
| Matched bin algebra | 2.512e-13 |

Maximum matched closure is 2.220e-16; maximum C2 closure is 1.137e-13.
Couplings and diagnostics are shared metadata, not quantities to sum over bins.

The next numerical work is a controlled lattice-error budget for endpoint bins.
Physical flavour/spin/boost content, external aligned comparisons and a
validated NNLO production-host adapter remain outstanding.


## Completed extended grid scan

The 16384-bin baseline scan (up to 65536 bins) completed successfully. The
frozen-source independent audit is in
`results/lepton_spectrum_grid_refinement_checked.json`, with executable, source,
configuration and raw-output hashes. All ten points pass component/bin closure,
continuum O(alpha) invariance, domain/compression and scale-input checks.
The additional shifts are dimensionless integrals:

| Change | Pair total | Largest pT-bin shift | Largest eta-bin shift |
| --- | ---: | ---: | ---: |
| Double FF panels | -1.002e-08 | 2.310e-08 | 8.612e-09 |
| Grid x2 | +2.578e-07 | 1.114e-05 | 2.248e-07 |
| Further grid x2 | +1.454e-07 | 4.462e-06 | 1.328e-07 |
| Double log span | -3.331e-16 | 4.885e-15 | 3.331e-16 |
| Tighten hard slivers | +1.271e-09 | 3.478e-09 | 8.796e-10 |

These are measured refinement shifts, not certified residual bounds. The
pT edge remains more sensitive than the pair total. The four correlated scale
variations also completed and are stored per bin in the checked JSON. This
remains the restricted rest-frame vector reference, separate from the new
physical-current NLO QCD decay pilots.
