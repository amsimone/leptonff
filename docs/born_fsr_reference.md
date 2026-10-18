# Independent QED FSR comparison on Born QCD

This benchmark compares bare muons using ReneSANCe/SANC fixed-order photonic
FSR, PHOTOS, PYTHIA, and leptonff at order alpha and matched LL/NLL. Production
is **Born QCD** throughout. This isolates the decay comparison and is distinct
from the DYTurbo NLO+NLL production campaign.

## Definitions and inputs

All tools receive the same SANC Born production distribution: pp at 13 TeV,
NNPDF31_nnlo_as_0118 member zero, factorization scale equal to the parent Q,
fixed alpha = 1/137.035999084, and muon mass 0.1056583755 GeV. The electroweak
Born inputs and widths are the pinned ReneSANCe defaults and are recorded in
each result's `settings`. No weak corrections, QED ISR, initial-final
interference, pairs, running alpha, QCD showers, or hadronization enter.
The common NNLO PDF is a benchmark input, not a claim of NNLO production.

* Z/gamma*: parent Q = 50–200 GeV; observed bare pair mass 50–116 GeV;
  no lepton pT or eta cuts. Visible mass bins are 0.25 GeV.
* W+ and W-: parent Q = 60–200 GeV; observed muon pT > 25 GeV and
  |eta| < 2.4, without neutrino-pT or transverse-mass cuts. Bins are
  0.25 GeV for muon pT, 0.05 for |eta|, and 0.5 GeV for transverse mass.
  The latter uses the final neutrino transverse momentum as particle-level
  missing momentum.

These are explicitly conditional parent domains for a validation study.
They are not a selected experimental measurement or a full-phase-space
phenomenology prediction. The parent upper bound exceeds the observed Z
window, so feed-in from 116–200 GeV is included. Feed-in above 200 GeV still
requires a domain study for a physical mass-window prediction. The plotted
finite pT/mT ranges do not exhaust the fiducial rate, which is accumulated
separately.

Each correction is (FSR-on / its paired Born) - 1. Ratios retain the full
within-bin channel covariance; independent runs are combined before division.
Ratios with a denominator below five standard deviations are withheld.
Shaded bands are one-standard-deviation **Monte Carlo errors only**.
Quadrature curves have refinement errors, not zero numerical uncertainty.

## Independent fixed-order reference

The optional ROOT-free driver compiles the public ReneSANCe 1.3.11 SANC
Fortran amplitudes and reuses the upstream hadron-process implementation.
The physics switches are iqed=4, iew=0, iqcd=0, with higher-order EW and
QED evolution disabled. LoopTools is initialized explicitly. The driver
respects SANC's different daughter ordering: W+ has the neutrino as leg 3,
whereas W- has the charged lepton as leg 3.

For W, massive real emission and virtual-plus-soft contributions are integrated
directly. Logarithmic energy and mass-regulated angular maps resolve soft and
collinear emission. Signed PDF weights are retained. At each parent Q an
independent integral of SANC's angle-integrated real spectrum, together with
its virtual/soft result, determines the inclusive correction. That correction
is subtracted times the Born distribution. This implements the same unitary
branching-fraction convention as leptonff/PHOTOS. The raw fixed-order result
can be recovered by adding the `inclusive_counterterm` channel.

For the mass-only Z measurement, a faster deterministic route integrates the
published SANC angle-integrated real expression in a whole-plus action:

    delta[F] = integral dGamma_real * (F(real) - F(Born)) / Gamma_Born.

The endpoint subtraction supplies the virtual contribution in the unitary
convention. This gives the complete normalized order-alpha mass distribution
without soft slicing. It uses SANC's massive real matrix element, not the
leptonff splitting or matching kernels. Born luminosities and the parent-mass
convolution use quadrature too. Muon masses remain finite in the reference;
the library's matching is at leading power.

There is an important upstream distinction: the differential NC real routine
obeys iqed, but its integrated `_1spr` counterpart contains ISR, FSR and
interference without consulting that switch. The pinned builder therefore
selects all eight FSR terms multiplying the current at s (`v0s`, `vas`),
preserving their coefficients. It excludes the currents at s' and the
interference terms. Generated GPL sources stay in the ignored build directory.
The direct differential MC calculation cross-checks this reduction.

The normalization convention and separation follow the
[SANC–PHOTOS study](https://doi.org/10.1140/epjc/s10052-013-2625-1).
ReneSANCe process controls are described in its
[hadron-mode paper](https://arxiv.org/abs/2207.04332).

## Shower and fragmentation configurations

PHOTOS 3.56 uses exponentiation and interference, no pairs, and the same fixed
alpha. Its full production-dependent ME correction is **off**, because this
adapter supplies a standalone decay record. The legacy W decay correction is
on. This is the existing universal-decay configuration, not a certification
of PHOTOS with a complete production record. The exponentiation tolerance is
1e-7 and the relative infrared cutoff is 1e-5.

PYTHIA 8.318 uses a forced timelike decay shower, fixed alpha, no photon
conversions or other showers, and decay matrix-element corrections on.
Individual daughter shower scales are set explicitly. The reader follows
`iBotCopyId()` through photon emissions and verifies that the measured muon
is final. A plain copy-chain lookup would incorrectly return the pre-FSR
lepton. Momentum conservation is checked in both shower interfaces.
The lepton shower cutoff is 1e-5 GeV; this has different units from PHOTOS'
relative cutoff. `--pythia-mec false` makes the separate uncorrected-shower
control; these settings must not be averaged together.

The library curves use the original order-alpha match and matched photonic
LL/NLL, with one subleading insertion globally, no pairs, and fixed alpha.
The W angular remainder uses common radiation seeds between accuracies.
The inclusive Z reduction now supports both LL and NLL and is checked
against the public two-leg operator for each choice.

## Two different matching comparisons

The primary output uses the library's original finite order-alpha matching.
It is the independent test of our first-order calculation. Fine W bins still
have appreciable angular-remainder Monte Carlo noise.

A separately labelled `--match-to-reference` output instead constructs

    SANC_NLO + (leptonff_resummed - leptonff_order_alpha).

The Born production and branching conventions are common. This is additive
matching to the independent massive NLO reference, with leading-power LL/NLL
corrections beyond NLO. It avoids sampling the noisy library first-order
term twice and includes the reference's finite-mass NLO terms. It is **a
different matching prediction**, not evidence that our original first-order
calculation agrees more precisely. The subtraction uses paired covariance;
the independent SANC error is then added. LL/NLL curves share that reference
error, so their band overlap is not a significance test.

## Reproduce

Install leptonff and provide an LHAPDF installation/PDF grid and a PHOTOS
source/library pair with the same C++ ABI. PYTHIA's `pythia8-config` and
gfortran must be on PATH. ROOT and Cuba are not required. LoopTools can be
built automatically from the pinned ReneSANCe archive or supplied with
`--looptools`.

```sh
python3 adapters/renesance/build.py \
  --lhapdf /path/to/lhapdf-prefix --leptonff /path/to/leptonff-prefix \
  --photos-source /path/to/photos/source --photos-library /path/to/libphotos.a

python3 adapters/renesance/run.py \
  --client build-renesance/leptonff_fsr_reference --pdf-data /path/to/PDFs \
  --output results/my-born-fsr --trials 2000000 --analytic-trials 100000 \
  --seeds 621 622 623 624 --hard-samples 4 --jobs 4

python3 adapters/renesance/compare.py results/my-born-fsr \
  --output results/my-born-fsr-plots
python3 adapters/renesance/compare.py results/my-born-fsr \
  --output results/my-born-fsr-shared-nlo --match-to-reference
```

The builder accepts `--archive` for an existing checksum-verified release.
It verifies extracted inputs against the archive before use and records source
and dependency hashes. Each campaign snapshots its executable so rebuilding
cannot change later seeds. Outputs include settings, commands, seed counts,
channel covariance, PNG/PDF/SVG plots, and spectrum TSV tables. Numerical
controls are kept separate by the analysis guards.
New campaigns namespace the input seed labels across processes and integration
modes. Additive matching of separately integrated MC results requires disjoint
streams; a single `all` run instead retains their joint covariance. For the
initial local campaign, the shared-NLO plot uses the disjoint 8-million-trial
W reference subset, while the original comparison combines all 16 million.

Refine Z with `--mass-nodes 32 --rapidity-nodes 128 --real-nodes 16
--ff-bins 1048576`. Use `--z-monte-carlo --modes fo` as a direct integration
control. Refine W soft slicing with `--omega 1e-6`; vary shower cutoffs and
MECs separately. Increasing W analytic statistics/hard samples is still needed
for precision original-matching spectra.

## Current local validation

The local campaign has 16 million fixed-order trials per W charge, 8 million
decays per shower/process, and 400,000 library production points with four
hard draws per accuracy and W charge. Z uses deterministic Born/FSR
quadrature, plus an 8-million-trial differential-MC cross-check.

* Refined SANC Z mass and library order-alpha curves differ by at most
  1.16e-6 of the Born bin rate. The coarse-to-refined quadrature shift is
  larger: 1.41e-5 for fixed order and 1.15e-5 for NLL. This is a leading-power
  validation, not a claim of one-part-per-million numerical certification.
* Direct Z MC versus the independent integrated expression has bin-pull RMS
  1.013. W soft-cutoff controls have pull RMS 1.101 (W+) and 0.987 (W-).
  These correlated-bin diagnostics are not global chi-square probabilities.
* Separate shower-cutoff runs (four million decays per process) reduce both
  technical cutoffs tenfold. Their bin-pull RMS values are 0.936–0.987 for
  PHOTOS and 1.022–1.034 for PYTHIA. No resolved cutoff bias appears at these
  MC uncertainties.
* Four-million-decay PYTHIA MEC-off controls are complete. For W+ the fiducial
  correction is -2.640 +/- 0.036 percent without MEC and
  -2.510 +/- 0.024 percent with MEC. For Z it is respectively
  -1.154 +/- 0.017 and -1.067 +/- 0.013 percent. These are MC errors;
  this comparison does not assign a theoretical uncertainty.
* Core mass-response/production/NLL tests and Python covariance/content guards
  pass. The final build also exercises automatic LoopTools construction.

The diagnostics can be regenerated for this documented campaign layout with
`python3 adapters/renesance/validate.py --root results/born-fsr`. The small
tracked snapshot is `benchmarks/born_fsr_reference.json`; its stochastic
values are not bitwise regression targets.

Original-matching results are under `results/born-fsr/comparison/`;
the separately matched results are under
`results/born-fsr/comparison-shared-nlo/`. The validation report is
`results/born-fsr/validation.json`. Generated data are intentionally untracked.

Remaining work: reduce fine-bin shower and W angular-remainder statistics;
test production-dependent PHOTOS ME corrections;
converge the parent domain for a selected fiducial measurement; and repeat the
validated decay comparison on resummed QCD production. Dressed finite-cone
NLL and exclusive multiphoton angular observables still require additional
theory/interface work. These Born-QCD plots do not complete those tasks.

## Leptonff-only accuracy and fragmentation-scale comparison

The subsequent [complete-radiator validation study](z_radiator_validation.md)
adds independent coefficient tests, a KKMC reference and a named matching
variation. Its external residuals are not included in the scale-only figures.

`compare_leptonff.py` isolates the library's own order-alpha, matched LL and
matched fixed-alpha photonic NLL calculations. It does not replace the finite
first-order coefficient with SANC. The LL label includes the one-loop boundary
and the LL-dressed hard coefficient; it is not evolution from a bare endpoint
with all finite matching omitted. NLL adds the conditioned two-loop photonic
kernel with one subleading insertion globally, as described in
[the photonic NLL contract](photonic_nll.md).

The optional runner argument `--frag-scale-factors 0.5 1 2` sets
`mu_frag = factor * pole_mass`. Both evolution and hard matching use the
varied scale; the physical parent Q remains in the hard logarithm. The
on-shell alpha, physical muon mass, initial scale mu0=m_mu, PDFs and QCD
production are held fixed. This is a QED fragmentation-scale scan, not a
QCD scale scan or a variation of the photon/pair content.

For example, the finest Z scan in the current local study is:

```sh
python3 adapters/renesance/run.py \
  --client build-renesance-v2/leptonff_fsr_reference \
  --pdf-data /path/to/PDFs --output results/z-accuracy/fine \
  --processes z --modes analytic --seeds 801 --jobs 1 \
  --frag-scale-factors 0.5 1 2 --ff-bins 2097152 \
  --mass-nodes 64 --rapidity-nodes 256 --real-nodes 16
```

Repeat into `grid-control` with 1,048,576 FF cells at the same quadrature,
and into `quadrature-control` with the finest FF grid and 32 mass / 128
rapidity nodes. Each directory contains all three scales. Then:

```sh
python3 adapters/renesance/compare_leptonff.py \
  --inputs results/z-accuracy/fine \
  --grid-control results/z-accuracy/grid-control \
  --quadrature-control results/z-accuracy/quadrature-control \
  --output results/z-accuracy/plots
```

The analysis rejects mismatched physics, mixed-resolution scale envelopes,
repeated deterministic calculations at a scale, and failed first-order
scale cancellation. Its ratio helper retains numerator/denominator covariance;
these Z inputs have no MC covariance because their integration is deterministic.
The compiled mass-response test additionally checks varied-scale first order
and agreement with the public two-leg matching operator at M_Z/2 and 2*M_Z.

Full-range and peak-zoom figures show the spectra, FSR/Born corrections,
resummed/order-alpha corrections and NLL/LL differences. Colored shading is
the three-point QED scale envelope. Hatching in the final panel separately
shows the sum of absolute grid and quadrature refinement shifts, maximised
over the three scales. This is an empirical numerical estimate, not a rigorous
error bound or a statistical band. The JSON/TSV also give those shifts for
every displayed ratio. Scale envelopes are not confidence intervals and
need not shrink monotonically or bracket the next logarithmic order.

The local nine-run study and figures are under
`results/born-fsr/z-accuracy-scales/`. W angular integration and shower
statistics are separate subsequent tasks.

The [small validation snapshot](../benchmarks/z_accuracy_scales.json) records
the input hashes, peak integrals and refinement checks. At the central scale,
LL changes the order-alpha rate in 90–92 GeV by +0.8634%; NLL changes LL by
+0.01381%. Differential resummation effects are larger than this integrated
peak average: LL/order-alpha spans roughly -3.28% to +2.17% across 50–116 GeV.
The NLL/LL change spans -0.0441% to +0.0250%. The numerical NLL/LL refinement
estimate is below one tenth of the increment in all 264 mass bins. These are
Born-QCD, unitary photonic benchmark results, not full-QED theory errors.
