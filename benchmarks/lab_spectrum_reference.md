# NLO W/Z observed-spectrum benchmark

Version 0.20 supplies the differential milestone described in the
[measurement/statistical contract](../docs/lab_spectra.md). It supports
observed W pT/absolute-eta/mT and charge-symmetric Z pair-mass/minimum-pT/
maximum-absolute-eta spectra on the MCFM NLO QCD production sample.

## Campaigns and validation

The nominal physics setup remains 13 TeV, W pole production and off-shell
Z/gamma*, pT(mu)>25 GeV, |eta(mu)|<2.4 and observed Z mass 60–120 GeV.
The available NNPDF31 NNLO set is retained for interface validation; the
matrix elements are NLO QCD. Coupled evolution has one unit-charge family
with its explicit MSbar input. PHOTOS uses fixed alpha and no pairs.

The coupled campaign contains 39 runs: three processes at central and zero
coupling, 8192 fraction bins, half/double fragmentation scale, smaller
geometric exclusions, independent seed 2089 and six noncentral QCD scale
points. Its explicit call multiplier is 0.05, so it is a validation/sensitivity
campaign rather than a precision production sample.

The aligned photonic campaign adds 12 runs: central, zero coupling, 8192 bins
and refined geometric exclusions for each process. Both tools use fixed alpha
and no pairs; the analytic result has LL evolution and one finite insertion
globally. This aligned comparison is distinct from the coupled NLL diagnostic.

The completed audits pass **54/54 coupled** and **24/24 photonic** required
checks. The coupled audit retains nine high-chi-squared host warning flags.
The photonic audit retains four corresponding W+ warning flags.
Required checks test implementation identities and numerical controls; their
success does not override those convergence diagnostics. Per-bin scale
envelopes, grid/cutoff/seed shifts and their sampling errors remain separate.

The final release executable additionally reproduces every existing rate
histogram and all 136 differential histograms bit for bit for W+, W− and Z
at the reduced central call budget. The generic batched geometry supports
both coupled full/first-order measurements and photonic evolved-minus-bare
measurements without changing the coupled reference.

All **55/55 core tests** pass. The suite adds literal differential tensor checks, shared-current
equivalence and an independent signed-event covariance/migration fixture.
The GCC host compiler also passes the new spectrum test. Active differential
bridge smoke checks cover 60 hard configurations per process in each branch,
including repeatability, finiteness, all eight component identities and every
bin-to-rate sum. The ordinary external PHOTOS conservation/radiation smoke
continues to cover 1000 points per process.

## Larger coupled central rates

These common-iteration rates use the full single-run call budget. Errors are
from grouped-event products, not subtraction of native histogram errors.
The new totals reproduce the previous one-family reference.

| Process | Matched [pb] | Matched minus continuum order1 [pb] | PHOTOS minus matched [pb] |
| --- | ---: | ---: | ---: |
| W+ | 5383.258 ± 22.881 | -0.0834 ± 0.7713 | -2.738 ± 6.313 |
| W− | 4126.390 ± 36.484 | 0.5498 ± 0.5866 | -37.964 ± 32.233 |
| Z/gamma* | 725.235 ± 1.813 | 0.0298 ± 0.0805 | -2.148 ± 0.982 |

The coupled PHOTOS differences above compare different coupling/pair models.
The aligned photonic central sample is recorded separately below; neither
comparison has physical multi-flavour NLL accuracy. Full per-bin numbers,
normalised-shape covariance matrices and plots are in each result root.

## Larger aligned photonic central rates

Both tools have fixed alpha and no pairs. The analytic matching retains one
finite insertion globally on LL evolution. Errors below are integration
statistics; the NLL physical-flavour and precision gates remain open.

| Process | Analytic NLO+LL [pb] | PHOTOS minus analytic [pb] |
| --- | ---: | ---: |
| wplus | 5382.983 ± 22.916 | -2.464 ± 6.398 |
| wminus | 4125.858 ± 36.598 | -37.432 ± 32.328 |
| z | 725.145 ± 1.817 | -2.058 ± 0.987 |

## Result locations

| Artifact | Location |
| --- | --- |
| Coupled 39-run audit and uncertainty components | `results/spectrum_controls_v020/report.json` |
| Aligned photonic 12-run audit | `results/photonic_spectrum_controls_v020/report.json` |
| Larger coupled central runs and spectrum reports | `results/spectra_v020/` |
| Larger aligned photonic central runs and reports | `results/photonic_spectra_v020/` |
| Release compatibility checks | `results/spectrum_release_regression_v020/report.json` |
| Final core tests | `results/spectra_release_ctest.log` |
| GCC spectrum test | `results/spectra_gcc.log` |
| Coupled/photonic direct-bridge smoke checks | `results/spectra_smoke.log`, `results/photonic_spectra_smoke.log` |

Each run contains `iteration_rates.csv` and `event_covariance.csv`. The latter
stores the upper-triangular event covariance for all 144 selected channels at
each production iteration. Analysis reports contain the full summed covariance,
normalised fractions, exact bin edges, native-histogram reconstruction checks,
chunk diagnostics and input-file hashes. Infinite overflow edges are `null`.
Plots are exported as PNG/SVG; signed bins and errors are retained.

## Reproduction

Build the isolated adapter as in the contract. Run the complete coupled
validation campaign in a new results root:

```sh
python3 adapters/mcfm/run_spectrum_controls.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-spectra/mcfm_leptonff --results results/new_spectra \
  --call-scale .05 --jobs 3
python3 adapters/mcfm/analyse_spectrum_controls.py results/new_spectra \
  --output results/new_spectra/report.json
python3 adapters/mcfm/run_spectrum_controls.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-spectra/mcfm_leptonff --results results/new_photonic \
  --model photonic --cases nominal zero grid8192 slivers --call-scale .05 --jobs 3
python3 adapters/mcfm/analyse_photonic_spectra.py results/new_photonic \
  --coupled-results results/new_spectra --output results/new_photonic/report.json
```

Use `--cases nominal --call-scale 1` in a separate root for the larger central
samples. Increasing calls does not replace grid/cutoff/seed studies. The full
physical-flavour/PDF/recoil precision gates are listed in the
[operational summary](../docs/operational_summary.md).

## Release provenance

The reproducible core and final adapter/analysis source snapshot is
`results/spectrum_source_v020/`, core digest
`ee82586a728660d35e5710e9c9a698dfd50ac44a77eaad00e49e2ef0949c6b7f`.
The final isolated executable is `build-mcfm-spectra-release/mcfm_leptonff`.
Per-run manifests retain the executable/source hashes actually used, including
the earlier coupled build; the release regression verifies bitwise agreement
of all rates, bins, raw iterations and event cross-products. Later analysis
changes are hashed separately and do not relabel the simulation source.
