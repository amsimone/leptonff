# Hadronic alpha and LO QCD reference (v0.22)

The [contract](../docs/hadronic_qcd_matching.md) states the approximation:
sourced OS hadronic input, explicit conversion to a fermionic MSbar coupling,
LO QCD evolution of quark/gluon parents and a pointlike light-quark boundary.
This is the optional combined-evolution study. Mixed kernels and fitted
hadronic fragmentation are not prerequisites to the separate pure-QED NLL
benchmark; the [accuracy audit](../docs/qed_accuracy_audit.md) corrects the
earlier overbroad gate language and defines the current comparison.

## Independent validation

The final AppleClang Release suite passes **60/60 tests**. Independent GCC
C++ and reference-solver checks pass as well. The new checks cover SU(3) LO per-parent momentum conservation, asymptotic
freedom, thresholded strong-coupling composition, the zero-alpha_s limit,
the HVP anchor without duplicate perturbative running and all-parent
QCD/QED lattice moments. QCD-dressed C1/C2 reproduce the actual FF and W/Z
laboratory match through a², leaving cubic remainders.

The independent Python check constructs LO QCD moments analytically and
uses ordered exponential actions with midpoint/Richardson refinement.
It uses independently continued QED source moments and independently evaluates
the hadronic conversion and leptonic running. It checks N=2/N=3 with QCD on
and off, all parent rows and the selected observed columns. The production
joint RK4 algorithm is not used to make the reference. Both compiler builds
agree with this reference within 6.1e-13; the all-parent lattice moment
difference is below 3.1e-14 and the anchor difference below 2.8e-17.
The external C/Fortran/PHOTOS ABI test checks 60 points per W+, W− and Z,
including all 136 spectrum columns and their rate/overlap identities.

The external Jegerlehner 2019 implementation is read from the SANC installation
and compiled without modifying it. Its SHA256 is
`b8b9694a8d5f2e218a51b911561b97f4b9a5f5e9c095c4d8fa9cb5ec7e2ec64f`.
It gives Delta alpha_had5,OS(MZ²)=0.027611194178462029 with combined error
0.00010873963335546854. The unpaired comparison to the PDG input is 0.076
standard deviations; this is a diagnostic between evaluations, not a fit
assuming independent hadronic datasets. Reports and the generated alternative
input are in `results/external_hvp_v022/`. The source routine's internal
interpolation is single precision.

## W±/Z controls

All **39 formal runs complete and all 101 required host checks pass**.
Thirteen reduced-call W+ runs retain the host high-chi-squared flag; W− and
Z have no such flags in this campaign. These flags remain in the report.

The formal campaign is `results/hadronic_controls_v022/`. There are thirteen
cases per process: nominal, zero alpha, QCD evolution off, hadronic input ±1
standard deviation, matching scales 1.5/3 GeV, grid 1024, alternative HVP input,
fragmentation scales M/2 and 2M, legacy one-family mode and the lepton-only mode.
Nominal cases use 512 bins, log span 64, signed conservative cells, four hard
samples, a common seed and the 0.05 call budget. Legacy regression retains its
4096-bin/span-40 setup. These are reduced-statistics and finite-grid controls.

The common-iteration nominal rates are:

| Process | Matched [pb] | Matched minus order 1 [pb] |
| --- | ---: | ---: |
| W+ | 5283.2048 ± 127.3097 | 0.77307 ± 0.25099 |
| W− | 4151.4665 ± 60.0169 | 0.50983 ± 0.24909 |
| Z | 715.4437 ± 8.4841 | -0.16929 ± 0.16129 |

These errors use grouped-event covariance and common iteration weights.
The differences use the paired higher-order column. Comparisons to v0.21
central values with a different grid/measurement response do not isolate a
hadronic effect; use the same-grid controls and their recorded numerical
sensitivities. There are no cross-run event products or significance claims
for the grid, input or model shifts. Native inverse-variance rates are a
different estimator and are not subtracted to obtain these numbers.

Every run is checked for spectrum coverage, rate/overlap sums, native-error
reconstruction and full grouped-event cross-bin/channel covariance. First-order
continuum/Born/PHOTOS invariance, zero-alpha identities with nonzero QCD,
matching-scale preservation of the high-scale anchor and input monotonicity
are audited. Legacy mode must reproduce all 33 rate histograms, 136 spectra,
raw iteration rates and event products bit for bit. The input ±1 standard
deviation controls cover only that input contribution.

Plots `{process}_nominal_spectra.png/.svg` label the pointlike boundary and
LO QCD explicitly. PHOTOS comparisons are diagnostic; its fixed coupling and
pair content differ from the analytic mode. Production uses the existing
NLO matrix elements with the NNLO PDF reference; the fragmentation strong
input is independently recorded. This is not a precision input/PDF fit.

## Reproduction

```sh
python3 python/check_external_hvp.py /path/to/hadr5x19.f \
  --output results/new_external_hvp
python3 adapters/mcfm/run_hadronic_controls.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-hadronic-final-v022-release/mcfm_leptonff \
  --alternative results/new_external_hvp/jegerlehner19.dat \
  --results results/new_hadronic --jobs 3
python3 adapters/mcfm/analyse_hadronic_controls.py results/new_hadronic \
  --legacy-reference results/flavour_matching_v021 \
  --output results/new_hadronic/report.json
```

Build the isolated adapter with the existing `build_adapter.py` recipe and
`--iteration-audit`, using the same C++ ABI as the host. Always use fresh
output roots. Input files and source/executable hashes are retained per run.


## Source and release provenance

The formal simulation used core source digest
`c5fdde13ed1ffe47f8203e3bec3fcdb75a9a1b03fba2672a8db6075771ddebd2`,
retained in `results/hadronic_simulation_source_v022/`. The final source is
`results/hadronic_source_v022/`, digest
`5c84cb0708bc82e291a8980195a6607bc4a03199d1070d88563275c0f77bde8a`.
The two subsequent changes correct a synthetic test geometry's zero-coupling
limit and reject calibration inputs starting above their matching scale.
The latter guard does not affect any valid campaign input. The exact patch is
`results/hadronic_release_delta_v022.patch`.

The final executable is `build-mcfm-hadronic-final-v022-release/mcfm_leptonff`.
Three fresh nominal W±/Z runs reproduce all 33 rates, 136 spectra, raw
iteration estimates and event cross-products bit for bit. The regression
report is `results/hadronic_final_regression_v022/report.json`. The legacy
mode separately reproduces the v0.21 reference bit for bit across the same
products. Actual per-run source/executable hashes are preserved without
relabeling earlier simulations.

The final validation summary and artifact index are
`results/hadronic_validation_v022.json` and
`results/hadronic_artifacts_v022.json`. Development probes and earlier failed
fixture-check logs are excluded from the final validation. The final core
suite and independent GCC checks use the corrected fixture.
