# Development history and reference commands

This preserves the detailed implementation milestones from the original README.
For the current entry point, see [README.md](README.md). Links to `results/`
refer to separately retained local validation artifacts, excluded from Git.

Version 0.25 packages the generator-independent decay/tensor response as an
installable C++ library and C/Fortran interface, with configurable bins,
grouped production archives and checksum-verified multi-model replay.
MCFM now uses the common spectrum evaluator and can export production
without evaluating QED. See the [package and replay guide](docs/packaging.md).
DYTurbo 1.4.2 native coupling is pending its source archive; the generic
spin/angular interface is ready, but no NNLO+NNLL result is claimed.
Package installation and W±/Z replay closure are recorded in the
[validation benchmark](benchmarks/package_reference.md).

Version 0.24 adds a separate **fixed-alpha photonic NLL** mode, with explicit
unitary conditioning of the published time-like photonic kernel and one
subleading insertion globally. W+/W-/Z adapter and comparison routing are
available as `--photonic-nll` / model `photonic_nll`. See the
[normalisation, matching and validation contract](docs/photonic_nll.md).
The frozen v0.23 high-level comparison plots still contain photonic LL;
this increment does not launch another production campaign.

The current target is **NLO QCD W/Z production with NLL final-state QED**,
plus a coupling/pair-aligned analytic-versus-PHOTOS comparison of W lepton pT
and |eta| and Z dilepton mass. Mixed QCD–QED fragmentation kernels and fitted
hadronic fragmentation are not prerequisites to this target. See the
[accuracy audit and comparison workflow](docs/qed_accuracy_audit.md).

The [completed comparison](results/qed_comparison_v023/comparison.md) provides
W± pT/|eta| and Z pair-mass plots, paired ratios, normalised shapes and
covariance exports. It passes 144 main consistency checks. W− uses a
ten-seed variance-control ensemble; the photonic domain-setting bug is fixed
with 48 bitwise default-regression checks. See the
[benchmark and remaining precision limits](benchmarks/qed_comparison_reference.md).

Version 0.22 adds sourced five-flavour hadronic alpha input, explicit MSbar
conversion/matching and LO QCD quark/gluon evolution. The new
`hadronic-pointlike` W±/Z benchmark keeps its light-quark boundary approximation
explicit. It is an optional combined-evolution study; missing mixed
coefficients and a fitted hadronic boundary do not block the pure-QED target.
It is not a full combined-accuracy data prediction.
All 60 core tests and 101 host checks pass across 39 control runs; thirteen
reduced-call W+ convergence flags remain recorded.
See the [contract](docs/hadronic_qcd_matching.md) and
[benchmark](benchmarks/hadronic_qcd_reference.md).

Version 0.21 adds full charged-flavour final-state QED evolution, finite
time-like threshold matching, leptonic on-shell/MSbar input conversion and
consistent W±/Z hard coupling and C1/C2 coefficients. The e/mu/tau preset
and strict custom charge/colour inputs are available alongside the previous
one-family and photonic modes. All 58 core tests and 90 required host checks
pass across 36 runs; twelve reduced-call W+ convergence warnings are retained. See the [flavour contract](docs/flavour_threshold_matching.md)
and [validation benchmark](benchmarks/flavour_threshold_reference.md).
Hadronic vacuum polarisation, QCD evolution of coloured parents, data-definition
alignment and a controlled differential precision budget remain open.

Version 0.20 completes observed-momentum W±/Z spectra on the NLO QCD host,
with full grouped-event bin/channel covariance, normalised shapes and PNG/SVG
plots. Both the one-family coupled NLL branch and the fixed-alpha photonic
PHOTOS benchmark are operational. The differential audit passes 54 coupled
and 24 photonic controls across 51 runs; host convergence flags remain visible.
See the [measurement contract](docs/lab_spectra.md),
[benchmark](benchmarks/lab_spectrum_reference.md) and
[operational summary](docs/operational_summary.md) for commands and remaining
physical-flavour, PDF, recoil and numerical/statistical precision gates.

Version 0.19 adds conservative signed integration over resolved fraction
cells for the coupled W/Z match, with the exact endpoint and consistent
C1/C2 overlap. It reduces artificial cut jumps and improves the measured
variance and grid response; the default atomic reference remains available.
See the [contract](docs/conservative_fraction_cells.md) and [benchmark](benchmarks/fraction_cell_reference.md). Physical flavours,
differential hard spectra and a complete uncertainty budget remain open.
The completed increment passes 53 core tests and 36 host checks; all 27
host runs completed. Gains vary across the two independent seeds.


Version 0.18 joins the spin-correlated W± and charge-symmetric Z hard
operators to coupled one-family MSbar NLL evolution, including photon-parent
conversion, charge-flipped pairs and the common running hard coupling.
The continuum NLO overlap and explicit second-order expansion are retained.
See the [contract](docs/coupled_lab_matching.md) and [benchmark](benchmarks/coupled_lab_reference.md). Physical flavours,
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
runs. See the [operator contract](docs/photonic_lab_dressing.md) and
[benchmark](benchmarks/photonic_dressing_reference.md). This is a fixed-alpha photonic finite-dressing
increment; full physical NLL still needs the coupled hard-parent/flavour/pair
extension and a controlled numerical/statistical budget.


Version 0.16 adds spin-correlated, unitary O(alpha) W and charge-symmetric
Z/gamma* decay matching on the NLO QCD host, with a fixed-alpha photonic
NLO+LL comparison to PHOTOS. See the [matching contract](docs/spin_decay_matching.md)
and [benchmark](benchmarks/spin_decay_reference.md). H1 is additive here;
its resummed dressing and physical multi-flavour NLL phenomenology remain open.

Version 0.15 supplies an operational **NLO QCD W±/Z host pilot with bare muons**:
laboratory FF convolution, an isolated MCFM adapter and a common-production
PHOTOS comparison. See the [contract](docs/lab_host_matching.md) and
[benchmark](benchmarks/nlo_host_reference.md). This is Born-parent FF-only
matching; version 0.16 extends it with the additive finite decay term.
The [W-first milestone](docs/w_nlo_milestone.md) now includes NLO Z validation
before a later NNLO QCD extension. Earlier standalone references are below.

Numerical references for [plan.md](plan.md):
distribution-aware perturbative lepton fragmentation, one-family MS̄ O(α)
boundary conditions, fixed-coupling photonic/non-singlet LL evolution, a
massless isotropic two-body decay at rest, and an angular-inclusive vector-current
hard+FF matching benchmark. The fiducial driver dresses the hard
coefficient with LL evolution on both muons and preserves full O(alpha) matching.
A separate NLO reference supplies the complete one-family lepton/photon matrix,
including pairs, and coupled real-moment evolution with two-loop MS̄ running.
Version 0.9 adds a signed, energy-weighted coupled x-space grid, endpoint and
domain controls, and a first-order-matched Born FF fiducial reference.
Version 0.10 supplies the subtracted photon-primary hard coefficient and its
fixed-fraction fiducial action. All hard-parent channels are joined with coupled
NLO evolution and a consistent hard coupling in an angular-inclusive pair-mass
moment reference. Version 0.11 completes the coupled fiducial hard convolution
for the same restricted theory, with all required hard-parent channels,
complete O(alpha) overlap restoration and numerical refinement controls.
See the [contract](docs/photon_pair_matching.md) and
[benchmark](benchmarks/pair_matching_reference.md).
The [coupled fiducial contract](docs/coupled_fiducial_matching.md) and
[benchmark](benchmarks/coupled_decay_reference.md) describe the new rate.
Version 0.12 adds the matched pair-mass spectrum, conservative mass-cell
integration, component-wise bin closure and correlated fragmentation/starting
scale scans. See the [spectrum contract](docs/mass_spectrum.md) and
[benchmark](benchmarks/mass_spectrum_reference.md).
Version 0.13 adds single-muon pT and signed eta spectra with the same pair cuts,
two entries per pair, full hard matching and correlated bin studies. See the
[contract](docs/lepton_spectra.md) and [benchmark](benchmarks/lepton_spectrum_reference.md).
Version 0.14 evaluates endpoint/resolved geometry directly on the FF lattice
and aligns a common cubic panel partition with pT/eta contacts. See the
[numerical contract](docs/endpoint_quadrature.md) and
[current benchmark](benchmarks/endpoint_quadrature_reference.md).

The fiducial LL reference retains one finite boundary or hard insertion globally.
It matches the complete O(alpha) overlap; its ingredients and accuracy limits
are explicit in the output.

The fixed-order and LL drivers output dimensionless acceptances. A separate
mass-only driver includes the finite O(α) hard coefficient and an exact massive
vector reference. The fiducial vector driver includes the pT/eta cuts and adds
the finite hard correction to the LL result. A further driver convolves that
hard correction with both LL-evolved muons. The coupled driver supplies the
one-family fiducial NLL match. Hard-dressed physical multi-flavour NLL hadronic
predictions remain outstanding. These are reference
calculations rather than matched NLO EW or NNLO Drell–Yan predictions.

## Build and check

Requires CMake ≥ 3.16 and a C++17 compiler; no external physics library is needed.
Python ≥ 3.8 is optional and enables an independent numerical reference check.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/leptonff_decay configs/reference.conf
./build/leptonff_ll configs/reference.conf
./build/leptonff_match configs/reference.conf
./build/leptonff_hard configs/vector_mass.conf
./build/leptonff_fiducial configs/reference.conf
./build/leptonff_dressed configs/reference.conf
./build/leptonff_kernels
./build/leptonff_moments
./build/leptonff_xspace configs/reference.conf
./build/leptonff_photon_pair configs/reference.conf
./build/leptonff_pair_moments configs/reference.conf
./build/leptonff_coupled configs/reference.conf
./build/leptonff_mass_spectrum configs/reference.conf configs/mass_bins.conf
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf
```

The driver prints JSON with the input parameters, approximation, compiler,
integration settings, moment residuals and tolerance-refinement shift. An explicit
configuration must supply every key; unknown, duplicated or invalid keys fail.

To reproduce the reference tables with independent numerical checks and
save source/executable/configuration hashes:

```sh
python3 python/check_reference.py build/leptonff_decay configs/reference.conf \
  --output results/reference.json
python3 python/check_ll_reference.py build/leptonff_ll configs/reference.conf \
  --output results/ll_reference.json
python3 python/check_matching_reference.py build/leptonff_match configs/reference.conf \
  --output results/matching_reference.json
python3 python/check_vector_mass_reference.py build/leptonff_hard configs/vector_mass.conf \
  --output results/vector_mass_reference.json
python3 python/check_fiducial_reference.py build/leptonff_fiducial configs/reference.conf \
  --output results/fiducial_reference.json
python3 python/check_hard_dressed_reference.py build/leptonff_dressed configs/reference.conf \
  --output results/hard_dressed_reference.json
python3 python/check_nlo_kernel_reference.py build/leptonff_kernels \
  --output results/nlo_kernel_reference.json
python3 python/check_moment_reference.py build/leptonff_moments \
  --output results/moment_reference.json
python3 python/check_xspace_reference.py build/leptonff_xspace configs/reference.conf \
  --output results/xspace_reference.json
```

The fixed-order independent script covers the supplied reference geometry;
it deliberately rejects unsupported cut geometries. The C++ driver supports
other parent-at-rest cuts, with analytic angular acceptance and explicit z
breakpoints. For a custom file, copy the config and change its numerical values.

The LL driver reports a baseline, doubled grid resolution, then a halved soft
cutoff at the finer resolution. The cutoff is a numerical regulator and its
effect must be checked. For explicit controls:

```sh
./build/leptonff_ll configs/reference.conf 262144 0.0002 32
```

The optional arguments are the number of lattice points, unresolved cutoff in
`u=-log(z)`, and total u span. Counts must be powers of two; the cutoff must be at
least half a cell. The driver also checks the refined configuration. It rejects
a domain with excessive periodic-wrap probability instead of silently aliasing
radiation into the endpoint. Defaults are the tested reference, not settings
validated for arbitrary scales or sharp measurements.

The LL reference uses a delta boundary, fixed α and only the non-singlet kernel.
It omits finite O(α) boundary terms included in `leptonff_decay`. Comparing those
two rates directly therefore does not isolate higher orders. The LL output
includes its own consistently truncated LL O(α) rate for that comparison.

`leptonff_match` adds the evolved signed boundary, retains at most one boundary
insertion across the two legs, and subtracts the actual lattice O(α) overlap
before adding the continuum FF coefficient. It accepts the same grid options.
This enables a consistent FF-level comparison to the full O(α) FF reference.
It still lacks the hard decay and the NLO splitting kernel required for NLL.

`leptonff_hard` reports separate hard and FF terms, their scale cancellation,
a decreasing-mass scan, and an exact massive reference for a pure vector current.
It accepts only a mass window: `pt_min=0` and an inactive eta cut (use 100).
It rejects the angular cuts in `reference.conf`. Its rates are normalised to
massless vector Born. This first-order benchmark is separate from resummation;
see [derivation and source audit](docs/vector_mass_matching.md) and the
[benchmark](benchmarks/vector_mass_reference.md).

`leptonff_fiducial` adds finite hard radiation with symmetric pT, eta and
mass cuts for an unpolarised vector at rest. It reports a massive reference,
hard-regulator and mass/scale scans, and three LL lattice refinements. The
hard coefficient is included once at O(alpha), without LL dressing. Its
coupling expansion reproduces the full O(alpha) result; it is not NLL.
It accepts the same optional lattice controls as `leptonff_match`. See the
[fiducial matching contract](docs/fiducial_vector_matching.md) and
[benchmark](benchmarks/fiducial_reference.md).

`leptonff_dressed` convolves H1 with pure LL on both muons, while retaining
at most one finite boundary or hard insertion globally. It keeps the full
O(alpha) match and adds the hard-generated subleading logarithmic towers.
It is still not NLL: P1 evolution, running and mixing are absent. It accepts the
existing lattice controls plus an optional FF quadrature panel count:
`config bins LL_cutoff log_span FF_panels`. The angular calculation costs
more than additive matching and reports separate FF quadrature, lattice,
LL cutoff and hard-geometry refinements. The mass-window configuration
also works. See [hard dressing](docs/hard_dressing.md) and the
[benchmark](benchmarks/hard_dressed_reference.md).

`leptonff_kernels` audits the complete MS̄ NLO lepton/photon matrix. It outputs
source-separated pieces, same/opposite-charge and valence moments, and finite
lower-bound actions. These kernels include real/virtual pairs; they are not
inserted into the photonic LL decay drivers. The valence combination conserves
charge after analytically cancelling the singular singlet. Every parent N=2
momentum sum is checked. See the [source audit](docs/nlo_kernel_audit.md) and
[photon derivation](docs/nlo_photon_kernels.md).

`leptonff_moments` solves the coupled real-N evolution with fixed, one-loop and
two-loop MS̄ running, retaining the O(alpha) FF boundary. It reports propagator
and FF matrices, running integrals and explicit first/second expansion
coefficients. Its input coupling is explicitly MS̄ at mu_0. It is a separate
reference for N>1; it is now checked against the separate x-space grid.
The independent checker uses analytic running and Magnus matrix integration.
See the [solver contract](docs/moment_evolution.md) and
[benchmark](benchmarks/moment_reference.md).

`leptonff_xspace` evolves the full matrix on an energy-weighted log-loss grid,
retaining the endpoint distribution and signed weights. It checks moment
recovery, two grid refinements and a doubled span at fixed spacing. Its Born
two-direction FF pair density subtracts the actual O(alpha) lattice overlap
and restores the continuum FF coefficient. It includes no hard H1 and is not
a full fiducial NLL rate. This driver's config alpha is explicitly interpreted
as MS̄ alpha(mu_0). The independent reference uses direct causal convolution
and a power series, with no FFT or RK4. The optional final argument sets the
baseline bins (default 4096). See the [contract](docs/xspace_evolution.md) and
[benchmark](benchmarks/xspace_reference.md).

## Current components

- `Distribution`: regular, delta and whole-kernel plus terms; signed linear
  combinations, finite lower bounds and real Mellin moments.
- `apply_pair`: tensor action with conditional cut breakpoints and both endpoint
  subtractions. `apply_pair_order1` keeps exactly O(α), avoiding an accidental
  product of two truncated FFs at O(α²).
- One-family `{mu-, mu+, gamma}` O(α) FF boundaries and LO time-like kernels;
  explicit first-order scale transport with starting-scale cancellation.
- Complete one-family NLO time-like lepton/photon matrix and valence kernel,
  with explicit pair interference and independent analytic/numerical moments.
- Coupled real-moment NLO evolution, two-loop MS̄ coupling, boundary transport
  and first/second coefficient extraction, checked with an independent solver.
- Coupled energy-weighted x-space evolution, linear-z endpoint projection,
  weighted wrap/tail diagnostics and FF-only fiducial overlap matching.
- Born four-vectors, collinear rescaling and isotropic fiducial acceptance.
- Mass-window and symmetric fiducial hard+FF matching, spherical-band angular
  averaging, exact massive vector radiators and an independent Dirac-trace current.
- Additive fiducial NLO+LL matching and LL dressing of H1 with complete O(alpha)
  overlap removal, signed FFT evolution and a reduced two-leg FF quadrature.
- `LLNonSinglet`: finite-cutoff, log-loss lattice evolution; full-kernel
  exponentiation, endpoint weight, continuum moment checks and numerical
  diagnostics. Symmetric two-leg decay convolution uses a prefix-sum reduction.
- `LLBoundaryNonSinglet` and `decay_ff_ll_match`: signed boundary transport,
  bilinear two-leg cuts, coefficient extraction and complete O(α) FF overlap.

See [conventions](docs/conventions.md), [kernel sources](docs/kernel_sources.md)
and [implementation status](docs/status.md). [benchmarks/reference.md](benchmarks/reference.md)
records the fixed-order reference. [LL evolution](docs/ll_evolution.md),
[backend decision](docs/backend_decision.md) and [LL benchmark](benchmarks/ll_reference.md)
document the new solver and its limitations.
See [boundary matching](docs/boundary_matching.md), its
[benchmark](benchmarks/matching_reference.md), and the
[NLO kernel audit](docs/nlo_kernel_audit.md) for the kernel work still outstanding.
