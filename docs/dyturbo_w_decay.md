# Differential DYTurbo W decay coupling

The optional `leptonff_dyturbo_decay` client now applies the common decay
functional to DYTurbo 1.4.2 **NLO+NLL W production**, separately in the
resummed, counterterm and V+jet contributions. It evaluates QED off, finite
order-alpha matching, photonic NLL and PHOTOS on shared production points.
This is an implemented coupling with pilot validation; a high-precision W
shape comparison remains an integration task. It does not support NNLO or Z
fiducial decay coupling yet.

## Physics and integration contract

The QCD configuration is the native NLO+NLL production benchmark: 13 TeV,
NNPDF31_nnlo_as_0118 member zero, and explicit hard/counterterm/V+jet orders
equal to one. DYTurbo calls this the primed hard-coefficient convention.
Its QED evolution is disabled. The client requires full Born angular support,
one production-domain bin and separate differential production contributions.
It never divides a correction by a Born event weight.

At each native point the public `DecayCalculator` includes the collinear
fragmentation response and finite angular remainder. The signed dual spin
projector retains the production angular information. The same linear
operator acts on positive and negative QCD terms. This is factorized QCD
production with matched QED decay, not a complete mixed QCD/EW calculation.

Defaults are W+ parent mass 60–100 GeV, muon pT > 25 GeV, |eta| < 2.4, no
neutrino-pT or transverse-mass cuts, and 48 |eta| bins of width 0.05. Decay
alpha is fixed at 1/137.035999084, m_mu = 0.1056583755 GeV and mu_frag =
80.385 GeV. Photonic NLL has one subleading insertion globally, finite
order-alpha matching and unitary inclusive normalization. Pairs, running
alpha, a QCD shower and production photons are absent.

PHOTOS 3.56 runs on an on-shell massive two-body decay record, with fixed
alpha, exponentiation, interference and its legacy W correction enabled.
The production-dependent full W/Z ME weights are explicitly disabled: the
record does not carry the incoming parton flavours/spin. The wrapper checks
that PHOTOS' actual ME channel stays zero and checks four-momentum closure.
The analytic and PHOTOS massive Born baselines are retained separately.
Small mass effects and finite-kernel differences must not be mistaken for
Monte Carlo errors. PHOTOS is process-global and is used serially in each
client process; independent processes are safe.

The original MCFM pilot's standalone W record could not activate PHOTOS'
production ME channel even when its flag was set. That earlier flag label
does not demonstrate an active full ME correction.

The client calls native 6D RES, **8D CT**, and 7D V+jet routines. The 1.4.2
steering value `intDimCT=6` selects a C++ routine that consumes eight random
coordinates, including two PDF convolutions. This distinction is checked
explicitly. The initial 6-coordinate CT pilot is invalid and is excluded.

Four VEGAS iterations use N/10, 2N/10, 3N/10 and 4N/10 trials. Their
estimates are combined with deterministic trial-count weights. Per-bin
5-by-5 covariance retains on/off/model correlations. The last-iteration
weight sum and variance are checked against Cuba. The fixed-budget Cuba
failure flag means the zero requested tolerance was not reached, rather
than certifying convergence. Cross-bin covariance is not in the compact
client result; the public production replay can calculate it.

## Build and run

Use a PHOTOS library compiled with the same C++ ABI as DYTurbo/leptonff.
On macOS, a GCC/libstdc++ PHOTOS library cannot be mixed with Clang/libc++.
No PHOTOS sources are bundled or modified.

```sh
cmake -S adapters/dyturbo -B build-dyturbo-interface \
  -DCMAKE_BUILD_TYPE=Release \
  -DDYTURBO_ROOT=/path/to/dyturbo/prefix \
  -DCMAKE_PREFIX_PATH=/path/to/leptonff/prefix \
  -DLEPTONFF_DYTURBO_DECAY_CLIENT=ON \
  -DCUBA_ROOT=/path/to/Cuba-4.2.2 \
  -DPHOTOS_SOURCE=/path/to/photos/source \
  -DPHOTOS_LIBRARY=/path/to/libphotos.a
cmake --build build-dyturbo-interface --parallel 4
python3 adapters/dyturbo/run_decay.py \
  --executable build-dyturbo-interface/leptonff_dyturbo_decay \
  --pdf-data /path/to/PDF/grids --output results/w-coupled-pilot \
  --calls 10000 --hard-samples 8 --photos-samples 32
```

The runner saves the input card, backend log, backend/PDF/client checksums,
native means/covariances, a precision assessment, and `production.lff`.
`precision_met=false` is an explicit result, not a reason to hide bins or
omit their covariance. Independent seeds and grid/geometry/cutoff controls
are needed in addition to the reported Monte Carlo errors.

The production archive contains QCD-only full-angle points. Each adaptive
iteration is an additive stratum, with its trial-count fraction already in
the stored weights. This matches the public grouped-integration convention.
New fragmentation grids, hard samples, cuts or observables can be replayed
within the recorded production support without rerunning QCD. The archive
does not provide arbitrary PDF/scale reweighting or extend its mass domain.

## Validation and remaining precision work

The local zero-alpha check closes all five predictions exactly bin by bin.
Export/replay reproduces Born/first-order/matched means within 2e-12 pb and
Born variances within 3e-10 pb². The native fiducial quadrature rate is
5084.70 +/- 6.59 pb (backend error estimate); the 10,000-trial-per-component
coupled pilot gives 4882 +/- 224 pb before QED, consistent within one sigma.
Its matched total on/off ratio is approximately 0.9784 +/- 0.0015.

These pilot numbers are not a precision W shape result. Independent sampling
of large CT and V+jet contributions, angular-remainder sampling, and sparse
fine-bin Born populations produce much larger errors in individual bins.
The next performance improvement should combine the finite QCD contribution
with correlated sampling or a production-tensor integration, while retaining
the differential client as a reference. Production archives also allow much
larger decay sampling without repeating the expensive Mellin inversions.
No final fine-bin W plot is certified by this milestone.
