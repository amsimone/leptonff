# Flavour and threshold reference (v0.21)

The [implementation contract](../docs/flavour_threshold_matching.md) specifies
the charged-flavour kernels, finite time-like maps, leptonic alpha conversion,
selected observed columns and W±/Z hard matching. The final isolated host
executable is `build-mcfm-flavours-v021/mcfm_leptonff`. The formal campaign lives
in `results/flavour_matching_v021/`; development probes and earlier builds in
other `flavour_*` directories are excluded from the release audit.

## Core and independent checks

The AppleClang Release suite passes **58/58** CTest checks. GCC 16 independently
passes the new flavour and threshold/evolution C++ tests against the host core.
Checks include all-parent momentum sums, N=1 valence cancellation, one-family
reduction, unlike-flavour charge/colour factors, the analytic time-like
matching moment, nonzero heavy matching at its mass, forward/backward alpha
transport and the regular zero-coupling limit of the running logarithm.

The complete matched FF and actual laboratory hard+FF remainders scale
cubically after subtraction through a². The latter checks W+, W− and the
charge-symmetric Z operator with both atomic and conservative-cell measures.
Threshold-scale changes leave no retained a² ln(Q) term in any parent row.
The FFT moments check all parents, not just the selected 3×3 subblock.

`python/check_flavour_reference.py` forms source moments from independent
space-like expressions plus time-like continuation, analytic charge factors
and an analytic threshold Mellin coefficient. It transports alpha by a
bracketed implicit solution and evolves by exponential midpoint with Richardson
extrapolation, independently of the production RK4 solver. Twelve cases cover
leptons and synthetic coloured charges, three matching-scale choices and
one-/two-loop running. N=2/N=3 lattice moments are checked separately against
the continuum result. The maximum evolution difference is below 6.4e-14,
the source-kernel difference below 2.6e-11 and midpoint-refinement diagnostic
below 2.3e-10, for both compilers. Detailed results are
`results/flavour_independent_v021.json` and
`results/flavour_independent_gcc_v021.json`.

The external PHOTOS/adapter smoke test passes 1000 points per W+, W− and Z
process. The new variable-flavour C/Fortran ABI check passes 60 points per
process, checking all 136 differential columns for repeatability, finiteness,
bin/rate closure and paired-difference identities.

## NLO QCD host campaign

All **36 runs complete and all 90 required checks pass**. Twelve reduced-call
W+ runs carry the host high-chi-squared warning; all are retained in the report.
W− and Z have no host convergence flags in this campaign. The W+ warnings recur across the common-seed controls and must be retained
when assessing statistical stability.

The formal campaign has twelve cases for each of W+, W− and Z:

- e/mu/tau with OS input and nominal tau matching at its mass;
- tau thresholds at m/2 and 2m;
- zero coupling;
- e/mu/tau with the same numeric input interpreted as MSbar;
- legacy one-family MSbar mode;
- custom lepton file with OS input;
- 8192 bins versus the nominal 4096;
- fragmentation scale M/2 and 2M;
- hard soft/collinear slivers 1e-6 versus 1e-5;
- a synthetic coloured-fermion pure-QED test profile with MSbar input.

Each uses signed conservative cells, four hard angular samples, one thread,
the same NLO production seed and a call scale of 0.05. Every spectrum report
reconstructs native errors and checks total-rate, overlap and full grouped-event
covariance closure. File/preset equality and v0.20 legacy regression require
bitwise equality of all 33 rate histograms, 136 spectrum histograms, raw
iteration estimates and event cross-products. Threshold and grid controls
preserve continuum NLO/Born/PHOTOS columns; fragmentation-scale cancellation
is required within 1e-5 pb per raw bin/iteration.

The nominal common-iteration totals are:

| Process | Matched rate [pb] | Matched minus order 1 [pb] |
| --- | ---: | ---: |
| W+ | 5281.9612 ± 127.2164 | -0.47060 ± 0.92108 |
| W− | 4150.9950 ± 59.9910 | 0.03826 ± 0.61129 |
| Z | 715.2659 ± 8.4332 | -0.34709 ± 0.45926 |

These are reduced-call validation outputs, not replacement precision central
predictions for the v0.20 larger samples. Errors use common iteration weights
and full grouped-event covariance; higher-order errors come from the paired
difference column. The native histogram's inverse-variance total/error is a
different estimator and is not subtracted to obtain these differences.

For the nominal preset, alpha_MSbar(mu0)=0.007358102740118272.
At M_W=80.385 GeV it is 0.007560092649193642; at M_Z=91.1876 GeV it is
0.007564691687917954. These include leptonic vacuum polarisation only.
The associated B logarithms are 22.774245523979573 and 23.27861084921563.

Each case has a JSON spectrum/covariance report. The three nominal reports
also have PNG/SVG plots named `{process}_nominal_spectra`. The campaign report
records all checks, raw-data hashes, bin sensitivities and host integration
warnings. Scale, grid and cutoff shifts are diagnostics rather than certified
errors; there are no cross-run event products or significance claims for
differences between these runs. PHOTOS comparisons are diagnostic because the
pair content and coupling convention differ.

## Reproduction

Build the core with the host's C++ ABI, then relink read-only host objects:

```sh
cmake -S . -B build-host -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 -DBUILD_TESTING=OFF
python3 adapters/mcfm/build_adapter.py /path/to/MCFM-10.3.1/build-gcc \
  --output build-mcfm-flavours-v021 --photos-build /path/to/photos \
  --iteration-audit
python3 adapters/mcfm/run_flavour_controls.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-flavours-v021/mcfm_leptonff \
  --results results/new_flavours --jobs 3
python3 adapters/mcfm/analyse_flavour_controls.py results/new_flavours \
  --legacy-reference results/spectrum_controls_v020 \
  --output results/new_flavours/report.json
```

Use a fresh results directory. The host root passed to the runners differs
from its `build-gcc` directory passed to the linker. Configurations, custom
flavour files and actual executable/source hashes are retained per run.

This completes general pure-QED flavour and threshold matching, including a
physical lepton preset. Hadronic vacuum polarisation and QCD evolution of
coloured fragmentation parents remain open. Extra-pair/data definitions,
PDF/alpha_s uncertainties and numerical/recoil precision gates remain listed
in the [operational summary](../docs/operational_summary.md).


## Release provenance

The frozen source is `results/flavour_source_v021/`, with core SHA256
`2a014098fbfbb688400e84819d5482f970c715090232856f35ffe466b4cb99ae`.
This equals the source digest recorded by the final host build and every
formal run. `results/flavour_validation_v021.json` records compiler/core,
independent-reference, ABI and campaign validation. The artifact hash index is
`results/flavour_artifacts_v021.json`. Simulation manifests preserve the actual
adapter files at launch; later documentation/plotting changes are hashed in
the release snapshot and artifact index without relabelling simulation sources.
