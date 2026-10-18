# Implementation status

For the current package layout, interfaces and scientific contract, see the
[version 0.26 manual](manual/index.md). The operational history below records
the preceding physics and generator milestones.

Version 0.25 packages the generator-independent decay/tensor response as an
installable C++ library and C/Fortran interface, with configurable bins,
grouped production archives and checksum-verified multi-model replay.
MCFM now uses the common spectrum evaluator and can export production
without evaluating QED. See the [package and replay guide](packaging.md).
DYTurbo 1.4.2 native coupling is pending its source archive; the generic
spin/angular interface is ready, but no NNLO+NNLL result is claimed.

Version 0.24 adds a separate **fixed-alpha photonic NLL** mode, with explicit
unitary conditioning of the published time-like photonic kernel and one
subleading insertion globally. W+/W-/Z adapter and comparison routing are
available as `--photonic-nll` / model `photonic_nll`. See the
[normalisation, matching and validation contract](photonic_nll.md).
The frozen v0.23 high-level comparison plots still contain photonic LL;
this increment does not launch another production campaign.

The active target is NLO QCD production with final-state NLL QED in an
explicit flavour/input prescription. The [accuracy audit](qed_accuracy_audit.md)
corrects the earlier overbroad requirements: mixed strong/QED splitting and
fitted hadronic fragmentation are not prerequisites to the pure-QED W/Z
benchmark. The completed comparison covers W pT/|eta| and Z pair mass, with
independent seeds and grid/domain/cutoff controls: 144 main consistency
checks, 15 hadronic-input checks and 48 bitwise adapter regressions pass.
The main campaign has no host convergence warnings; five reduced-call W+
warnings remain in the auxiliary HVP-input study. W− uses ten independent
seeds at QCD dipole alpha=0.01 to reduce a measured production-variance
problem. The photonic adapter now honours its requested domain. See
[results](../results/qed_comparison_v023/comparison.md) and
[precision limits](../benchmarks/qed_comparison_reference.md).

Version 0.22 adds sourced five-flavour hadronic alpha input, explicit MSbar
conversion/matching and LO QCD quark/gluon evolution. The new
`hadronic-pointlike` W±/Z benchmark keeps its light-quark boundary approximation
explicit. Missing mixed coefficients and a fitted hadronic boundary apply to
that optional combined-evolution extension, not the pure-QED benchmark.
This is not a full combined-accuracy data prediction.
All 60 core tests and 101 host checks pass across 39 control runs; thirteen
reduced-call W+ convergence flags remain recorded.
See the [contract](hadronic_qcd_matching.md) and
[benchmark](../benchmarks/hadronic_qcd_reference.md).

Version 0.21 adds full charged-flavour final-state QED evolution, finite
time-like threshold matching, leptonic on-shell/MSbar input conversion and
consistent W±/Z hard coupling and C1/C2 coefficients. The e/mu/tau preset
and strict custom charge/colour inputs are available alongside the previous
one-family and photonic modes. All 58 core tests and 90 required host checks
pass across 36 runs; twelve reduced-call W+ convergence warnings are retained. See the [flavour contract](flavour_threshold_matching.md)
and [validation benchmark](../benchmarks/flavour_threshold_reference.md).
Hadronic vacuum polarisation, QCD evolution of coloured parents, data-definition
alignment and a controlled differential precision budget remain open.

Version 0.20 completes observed-momentum W±/Z spectra on the NLO QCD host,
with full grouped-event bin/channel covariance, normalised shapes and PNG/SVG
plots. Both the one-family coupled NLL branch and the fixed-alpha photonic
PHOTOS benchmark are operational. The differential audit passes 54 coupled
and 24 photonic controls across 51 runs; host convergence flags remain visible.
See the [measurement contract](./lab_spectra.md),
[benchmark](../benchmarks/lab_spectrum_reference.md) and
[operational summary](./operational_summary.md) for commands and remaining
physical-flavour, PDF, recoil and numerical/statistical precision gates.

Version 0.19 adds conservative signed integration over resolved fraction
cells for the coupled W/Z match, with the exact endpoint and consistent
C1/C2 overlap. It reduces artificial cut jumps and improves the measured
variance and grid response; the default atomic reference remains available.
See the [contract](conservative_fraction_cells.md) and [benchmark](../benchmarks/fraction_cell_reference.md). Physical flavours,
differential hard spectra and a complete uncertainty budget remain open.
The completed increment passes 53 core tests and 36 host checks; all 27
host runs completed. Gains vary across the two independent seeds.


Version 0.18 joins the spin-correlated W± and charge-symmetric Z hard
operators to coupled one-family MSbar NLL evolution, including photon-parent
conversion, charge-flipped pairs and the common running hard coupling.
The continuum NLO overlap and explicit second-order expansion are retained.
See the [contract](coupled_lab_matching.md) and [benchmark](../benchmarks/coupled_lab_reference.md). Physical flavours,
thresholds, differential hard spectra and a controlled phenomenology error
budget remain open. PHOTOS comparisons in this coupled branch are diagnostic
because the coupling scheme and pair content differ.
The completed audit passes 52 core tests, 486 independent polarised
RECOLA points and all 30 host controls. The small higher-order effects
remain limited by integration/grid resolution.


Version 0.17 dresses the spin-correlated W± and charge-symmetric Z hard
operator and the MSbar boundary with photonic LL, retaining one finite
insertion globally. The exact NLO overlap and previous rate columns are
preserved. Common-iteration rate/covariance audits accompany the NLO QCD
runs. All 51 core tests and 21 new host controls pass. See the [operator contract](photonic_lab_dressing.md) and
[benchmark](../benchmarks/photonic_dressing_reference.md). This is a fixed-alpha photonic finite-dressing
increment; full physical NLL still needs the coupled hard-parent/flavour/pair
extension and a controlled numerical/statistical budget.


Version 0.16 supplies leading-power spin-correlated W hard+FF matching and
the parity-even Z/gamma* counterpart for charge-symmetric rates. An optional
MCFM angular-tomography estimator preserves signed NLO QCD weights without
Born divisions. Fixed-alpha photonic NLO+LL and first-order/exponentiated
PHOTOS branches have explicit normalisation, pair and proposal conventions.
See [contract](spin_decay_matching.md) and [benchmark](../benchmarks/spin_decay_reference.md).
The W massive anchor, real-current spectra, Ward identities and spin limits
are tested; 50/50 core tests pass. Independent RECOLA tree checks cover 486
polarised points, agreeing within 3.2e-14. Full hard-dressed NLL, physical
multi-flavour content and a
phenomenology error budget remain open. The iteration audit exposes larger
common-weight errors for W- than the native per-histogram errors; the initial
projected-component statistical flag is retained alongside the passing
physical-rate controls.


Version 0.15 implements the first NLO QCD W±/Z hadronic FF-only pilot. A working
isolated MCFM 10.3.1 adapter applies one-/two-leg laboratory measurements before
observed cuts, retains signed real/dipole grouping and restores the FF O(alpha)
overlap. Native/zero-coupling checks, W±/Z grid and dipole-alpha controls and
same-production PHOTOS 3.56 comparisons for W+, W− and Z/gamma* are recorded in
the [contract](lab_host_matching.md) and [benchmark](../benchmarks/nlo_host_reference.md).
All 48 CTest checks pass. That earlier FF-only pilot does not provide a physical full-NLL prediction;
the subsequent additive hard match is described above. NLO W/Z work precedes the NNLO QCD extension.

The first executable batch implements plan steps 1–4 at reference level:

- C++17/CMake core with no third-party physics dependencies.
- Distribution action, finite bounds, moments and two-leg convolution.
- One-family MS̄ O(α) boundaries, LO time-like kernels and first-order scale
  transport; analytic moments, charge conjugation and momentum checks.
- One-family MS̄ NLO time-like charged-lepton block, source-separated pair
  contributions and a conserving valence combination, independently checked.
- Complete one-family NLO photon channels and per-parent momentum sums;
  independent space-like-plus-continuation checks of time-like formulas.
- Coupled real-moment NLO evolution with fixed, one- and two-loop coupling,
  O(α) boundaries, O(α²) coefficients and an independent Magnus reference.
- Coupled signed x-space lattice, energy tilt, endpoint projection, weighted
  wrap/tail controls and independently checked Born FF fiducial matching.
- Parent-at-rest Born decay, fiducial measurement and exactly O(α) FF action.
- Independent Python Simpson reference, provenance output and reproduction commands.
- Fixed-coupling photonic/non-singlet LL evolution, with the numerical endpoint
  atom, moment/expansion/composition checks and independent Poisson reference.
- Fast two-leg LL decay convolution and separate grid/soft-cutoff studies.
- LL dressing of H1 on both muons, signed hard-cell FFT transport, reduced
  moment-preserving FF quadrature and independent Poisson/angular checks.
- Symmetric fiducial vector hard+FF matching, exact rotational acceptance,
  massive reference and additive NLO+LL matching with expansion checks.
- Angular-inclusive vector-current finite hard coefficient plus FF matching,
  exact massive reference, independent Dirac trace and soft-slicing check.
- Evolved signed O(α) diagonal boundary, with at most one insertion across both
  muons; complete O(α) FF-only additive overlap matching and independent checks.

Stage A now has fixed-order and all-order non-singlet LL reference calculations,
with explicit numerical regulators and refinement controls. The LL backend is
standalone; [backend decision](backend_decision.md) records its scope. No external
evolution package has been forked. The coupled reference NLL backend is now
standalone; its physical-content and production-host extensions remain open.
No host or external massive-decay code access is assumed; no authors have been contacted.

Version 0.3 supplies the boundary-generated subleading terms and matches the
full O(α) FF coefficient. Version 0.4 adds hard+FF matching only for the separate
angular-inclusive vector mass spectrum. Version 0.5 extends that calculation
to symmetric fiducial cuts and additive NLO+LL. Version 0.6 dresses H1 with
LL on both muons. It is not full NLL: `P1`, running/mixing and physical
spin-correlated hard decay remain missing. [Boundary matching](boundary_matching.md) defines the prescription;
[NLO kernel audit](nlo_kernel_audit.md) records source equations and open mapping.

Version 0.7 implements the separate NLO lepton block, including opposite-charge
pair interference. Its charge/valence moment vanishes without altering the
source delta terms; the source's photonic component alone has a nonzero first
moment. The new `leptonff_kernels` audit driver has independent numerical and
analytic checks. The LL decay drivers still omit P1 evolution.

Version 0.8 completes the one-family NLO matrix and its real-moment running
evolution. The new `leptonff_moments` driver checks fixed, one-loop and two-loop
coupling, boundary transport and explicit O(α²) coefficients. Matrix composition,
ordering, backward inversion and momentum tests pass. The independent solver
uses analytic running and Magnus matrix steps. This is a moment reference;
endpoint-controlled x-space evolution, hard-coupling alignment and extra-pair
selection remain required before a fiducial NLL prediction. See the
[solver contract](moment_evolution.md) and [benchmark](../benchmarks/moment_reference.md).

Version 0.9 supplies the coupled x-space reference with linear-z endpoint
projection, exact source N=2/N=3 kernel moments up to tail omission, and higher
moment/grid convergence. The independent check uses causal power-series
convolution. A Born two-direction inclusive identified-pair FF action is matched
through O(alpha); it contains no hard H1. See the
[x-space contract](xspace_evolution.md) and [benchmark](../benchmarks/xspace_reference.md).

## Next work

Version 0.10 supplies the subtracted photon-primary hard coefficients, their
fixed-fraction fiducial action, and an all-hard-parent inclusive pair-mass
moment match. It joins coupled NLO evolution and alpha_MSbar(M) from one common
MSbar input, specifies inclusive distinct-parent pair counting, and verifies
cancellation of the retained alpha² mass-log scale dependence. See
[contract](photon_pair_matching.md) and
[benchmark](../benchmarks/pair_matching_reference.md).

Version 0.11 joins those ingredients for the coupled fiducial rate in the same
one-family unpolarised-vector model. All hard-parent channels and the common
MSbar coupling are included. Causal hard projections and a linear signed cubic
FF quadrature permit complete numerical overlap removal and explicit second-order
extraction. Independent causal/azimuth checks and grid, panel, sliver and scale
scans validate this restricted-model matched reference. See
[contract](coupled_fiducial_matching.md) and
[benchmark](../benchmarks/coupled_decay_reference.md). Physical spin/current,
multi-flavour and host gates remain open.

Version 0.12 adds the matched opposite-charge pair-mass spectrum. Shared
cumulative operators partition every component and expansion coefficient;
conservative resolved loss-cell integration controls mass-edge aliasing while
retaining the Born atom. Per-bin grid/FF/sliver refinements and correlated
fragmentation/starting-scale scans are recorded. Independent source/literal
bin and continuum O(alpha) checks pass. See [contract](mass_spectrum.md) and
[benchmark](../benchmarks/mass_spectrum_reference.md).

Version 0.13 adds both-muon pT and signed eta histograms with unchanged pair
selection. Mean one-leg cumulatives carry the full hard/FF match, common source
support, and first/second expansion; every bin sum is twice the pair rate.
Independent literal lattice and allowed-azimuth checks, selected massive
comparisons, and bin-specific numerical/scale studies are recorded in the
[contract](lepton_spectra.md) and [benchmark](../benchmarks/lepton_spectrum_reference.md).
Version 0.14 evaluates endpoint/resolved geometry directly on the FF lattice
and splits a common FF partition at the pT and eta contacts. This prevents
interpolation across forbidden pT support and removes panel interpolation from
the hard C2 coefficient. Signed lattice/support and expansion regressions and
a larger independent endpoint/contraction audit accompany the new refinements.
See the [numerical contract](endpoint_quadrature.md) and
[current benchmark](../benchmarks/endpoint_quadrature_reference.md).
At 8192 FF bins the largest four-to-eight-panel pT shift falls from 1.16e-4
to 3.25e-8; C2 has no panel change. The last FF-grid doubling still changes
a pT bin by 1.11e-5. These measured shifts do not constitute error bounds.
The broader physical uncertainty budget remains open.

1. Quantify the remaining FF-grid and resolved-interior errors after stabilising
   endpoint pT panel integration. Extend the recorded mass/pT/eta refinements
   and correlated scale studies into a fuller numerical and physics uncertainty
   budget. Unknown two-loop constant matching terms and missing physical content remain explicit limitations.
2. Define physical active flavours, thresholds and coupling inputs; maintain
   explicit inclusive pair counting. Exclusive two-muon selection is outside
   the current single-particle FF prescription.
3. Extend the coupled fiducial hard decay to the intended spin/current and
   boosted-parent content. The physical G2 gate remains open; see
   [hard dressing](hard_dressing.md).
4. Extend the operational NLO W/Z FF-only adapter with finite W hard matching,
   then Z hard matching, physical input alignment and observed pT/mass/mT bins.
   Complete differential refinement controls and full covariance.
   Compare against PHOTOS/PYTHIA using common production before promoting
   experimental-data or NNLO-host comparisons.
5. Establish external comparisons with aligned production, coupling and pair
   conventions before making a physical precision claim.

Fast C++ tests cover analytic distributions, a nonseparable cut with exact
double-plus action `-pi²/6`, moment identities, scale transport, kinematics,
acceptance and explicit O(α) truncation. The Python reference uses a different
quadrature and rationalised endpoint integrand to check the fiducial coefficient.
LL checks add continuum moments, direct composition, first/second expansion
coefficients, literal two-leg cut sums and an independent non-FFT Poisson solver.
Matching checks add signed boundary cells/convolution, a coupling scan proving
the O(α) overlap is removed, the no-evolution limit, global boundary truncation
and an independent signed Python convolution/matching calculation.

The vector benchmark checks finite hard moments, hard+FF scale cancellation,
starting-scale cancellation, convergence to the massive reference as m_mu/M
shrinks, and matrix-element gauge/kinematic identities. Independent Python
virtual-plus-soft slicing checks the inclusive normalisation. The reference
Release build uses AppleClang 21 and 53-bit Real.

Version 0.5 checks the spherical-band formula against azimuthal integration,
the exact massive scalar current against explicit Dirac traces, inclusive
reductions, eta-only and real-tail cuts, and massive/leading-power convergence.
It tests full fixed-order scale cancellation, the alpha-squared onset of the
matched remainder, no-evolution and zero-coupling identities, and configuration
compatibility. Independent Python angular quadrature checks the finite geometry
and complete O(alpha) coefficient. [Results](../benchmarks/fiducial_reference.md)
record numerical slivers and LL grid/cutoff convergence separately.

Version 0.6 adds signed hard-cell and literal convolution checks, the mass-only
U(2*t) reduction, polynomial FF quadrature moments, and alpha-squared onset of
the hard-dressing correction. Independent Python sums Poisson emissions and
signed convolutions without FFTs, then checks the dressed geometry using
azimuthal quadrature. [Results](../benchmarks/hard_dressed_reference.md) report
the larger FF quadrature refinement separately from the smaller grid/regulator
changes. Version 0.7 adds lepton-block moments, finite-bound actions, analytic
singlet cancellation and explicit rejection of unavailable photon channels.
Version 0.8 adds the complete matrix, running moment evolution and independent
continuation/Magnus checks. Version 0.9 adds weighted x-space projection,
independent causal evolution and fiducial FF matching. The version 0.9 reference
Release build passed all 31 checks. The [lepton-block benchmark](../benchmarks/nlo_kernel_reference.md) records
the 0.7 milestone; the [moment benchmark](../benchmarks/moment_reference.md)
records the moment checks; the [x-space benchmark](../benchmarks/xspace_reference.md)
records the grid/domain refinements. Version 0.10 adds photon-parent massive,
fiducial and scale checks and the all-hard-parent matched moments; its Release
build passed 37 tests. The [pair-matching benchmark](../benchmarks/pair_matching_reference.md)
records those checks. Version 0.11 adds coupled fiducial hard matching,
first/second expansion and mass-log scale checks, independent causal/geometry
contraction and separate numerical refinements. All 40 CTest checks pass in
the version 0.11 reference Release build. Version 0.12 adds mass-bin endpoint,
direct-tail, conservative cell, component/expansion closure and scale checks;
all 43 CTest checks passed at that milestone. Version 0.13 adds single-muon
mean-selection, signed non-endpoint tensor, Jacobian/contact, asymmetric eta,
two-entry closure, selected massive and per-bin scale checks. Independent
literal-bin and selected hard-geometry checks pass. All 46 CTest checks pass in
the version 0.13 Release build. Version 0.14 adds exact-contact ownership,
forbidden-support, literal mixed muon/photon tensor and endpoint-expansion
checks. Its independent audit includes direct endpoint angles and C2
reconstruction. All 47 checks pass in the current Release build.
