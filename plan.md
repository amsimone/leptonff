# Bare-muon W/Z production with lepton fragmentation — implementation plan

Proposed research extensions (18 October 2026): the
[extension notes](docs/manual/extensions.md) distinguish neutral-current NNLO
QED + collinear NLL matching from the dedicated W calculation, endpoint
resummation from the existing endpoint quadrature, and coherent soft/collinear
matching from cone-plus-veto dressed-electron predictions. They specify
overlap subtraction, physical veto efficiency, pair/axis definitions and
validation gates. These are proposals; no new resummation or production run
is enabled by the notes.

Architecture implementation (0.27): analytic projected hard actions, exact
fraction-pair stencils, reusable compiled measurements, deterministic threaded
replay, full/block/selected/diagonal covariance, shared C/Python measurement
views, explicit matched coupling inputs, batch C/Fortran/Python APIs, and
persisted tensor responses at explicit production nodes. See the
[performance and reuse guide](docs/manual/performance.md). Response interpolation
with an independent error budget and the native DYTurbo production adapter
remain separate work. These changes do not launch a new QCD campaign.

Package refactor (0.26): standalone fragmentation/evolution core, explicit
alpha inputs, distribution convolution/bin APIs, numerical Python bindings,
installed Fortran module, scientific manual and canonical benchmarks. See
[the current manual](docs/manual/index.md). The physics development history
and outstanding phenomenology goals below remain applicable.

Version 0.25 packages the generator-independent decay/tensor response as an
installable C++ library and C/Fortran interface, with configurable bins,
grouped production archives and checksum-verified multi-model replay.
MCFM now uses the common spectrum evaluator and can export production
without evaluating QED. See the [package and replay guide](docs/packaging.md).
DYTurbo 1.4.2 native coupling is pending its source archive; the generic
spin/angular interface is ready, but no NNLO+NNLL result is claimed.

Version 0.24 adds a separate **fixed-alpha photonic NLL** mode, with explicit
unitary conditioning of the published time-like photonic kernel and one
subleading insertion globally. W+/W-/Z adapter and comparison routing are
available as `--photonic-nll` / model `photonic_nll`. See the
[normalisation, matching and validation contract](docs/photonic_nll.md).
The frozen v0.23 high-level comparison plots still contain photonic LL;
this increment does not launch another production campaign.

Version 0.17 dresses the spin-correlated W± and charge-symmetric Z hard
operator and the MSbar boundary with photonic LL, retaining one finite
insertion globally. The exact NLO overlap and previous rate columns are
preserved. Common-iteration rate/covariance audits accompany the NLO QCD
runs. See the [operator contract](docs/photonic_lab_dressing.md) and
[benchmark](benchmarks/photonic_dressing_reference.md). This is a fixed-alpha photonic finite-dressing
increment; full physical NLL still needs the coupled hard-parent/flavour/pair
extension and a controlled numerical/statistical budget.


Implementation update, v0.16: finite spin-correlated W and parity-even,
charge-symmetric Z hard-decay operators now exist, with an additive MCFM NLO
QCD rate adapter and a fixed-alpha, pair-free photonic NLO+LL branch for
PHOTOS comparisons. Their [contract](docs/spin_decay_matching.md) separates
unitary decay normalisation, spin tomography, finite-remainder cutoffs and
first-order generator proposal sampling. Before phenomenology, close the
host controls and numerical/statistical error budget, dress H1 on all hard
parents for NLL, extend physical flavour/pair content, and add radiative
observable bins and widen the independent polarised matrix-element comparison.


17 October 2026 · Simone Amoroso · revised for implementation

## Goal and scope

Build a perturbative lepton-fragmentation library and an observable-level convolution layer. The first hadronic milestone is pp → W± → μ±ν at NLO QCD, with bare-muon fiducial rates and differential comparisons against PHOTOS/PYTHIA on common production. Start with an explicitly labelled FF-only integration pilot, then supply the W-specific O(α) hard decay and validate the matched prediction. Test the shared machinery for NLO QCD Z/gamma* cross sections next, before extending to NNLO QCD production and mμμ, single-muon pT and η spectra. Deliver reproducible numerical errors, scale studies, a documented library and at least one validated host adapter. Experimental-data comparisons follow the aligned tool comparisons.

### First hadronic milestone: NLO QCD W → μν

This milestone takes priority over NNLO host integration. Its contract and
initial local MCFM source probe are in [the W milestone](docs/w_nlo_milestone.md).
The existing neutral-current references remain validation fixtures.

Version 0.15 supplies the W± and Z/gamma* FF-only MCFM adapter, laboratory
actions and common-production PHOTOS comparisons. Native identity, zero-alpha,
W±/Z grid and subtraction-parameter controls are available. At that stage the finite hard
decay match was absent: see the [FF-only contract](docs/lab_host_matching.md)
and [benchmark](benchmarks/nlo_host_reference.md). Differential controls and
the physical-input/accuracy gates remain open.

1. Probe the locally available MCFM 10.3.1 W± processes at fixed-order NLO QCD
   with leptonic decay and QED corrections disabled. Reproduce its native
   rates/histograms with identity fragmentation before adding FFs. Preserve
   signed real/subtraction contributions, spin correlations and QCD recoil.
2. Implement the one-leg laboratory measurement with `p_mu = z*k_mu` and
   unchanged neutrino momentum. Define truth missing transverse momentum from
   the neutrino; visible collinear radiation must not become missing energy.
   Apply observed-muon, neutrino and transverse-mass cuts after fragmentation.
3. Compare pure NLO QCD, expanded O(α) FF action and resummed LL on the same
   production weights. This is an FF-only pilot, with its finite-hard and
   physical-content limitations explicit; it is not the final matched NLL result.
4. Supply a spin-correlated, subtracted W → μν(γ) O(α) hard decay and a massive
   or independent fixed-order check, including the gauge-consistent charged-W
   radiation contribution. Match its actual overlap to the FF evolution and
   define the width/branching normalisation and physical coupling/flavour inputs.
5. Compare bare-muon fiducial W+ and W− rates and eta spectra first, followed
   by coarse pT and transverse-mass spectra. Use an aligned independent
   calculation or FSR model; select experimental tables only after confirming
   their post-FSR bare definition. Record QCD/PDF, fragmentation and numerical
   uncertainties separately. Defer a W-mass extraction and narrow Jacobian bins.

The first host candidate is supported by local source/executable presence and
documented NLO W support, not yet by a completed leptonff adapter run. There is
one charged fragmentation leg, but the finite W decay cannot be obtained by
halving the neutral-current hard coefficient. The running neutral-current grid
study continues independently; its numerical conclusions are not W validation.

Per-mille control is a **target to establish**, not a consequence of the NLL label or a small scale band. Separate numerical precision, collinear logarithmic accuracy, fixed-order completeness, pole-approximation errors and omitted physics in every result. A publication is an outcome of the calculation; a visible disagreement with Photos is not an acceptance criterion.

At the planning baseline the workspace contained only this plan. Standalone
fixed-order, non-singlet LL, FF-boundary-matching and angular-inclusive vector
hard-matching references are implemented. A fiducial unpolarised-vector
hard+FF calculation, additive NLO+LL match and LL dressing of H1 on both
muons are now available; see
[implementation status](docs/status.md).
A separate complete one-family NLO lepton/photon matrix and coupled real-moment
evolution with two-loop MS̄ running are also available. Independent continuation,
moment, momentum and O(α²) checks pass. Version 0.9 adds the endpoint-controlled
coupled x-space grid and FF-only fiducial matching.
Version 0.10 supplies the hard photon-primary coefficients and a matched
inclusive pair-mass moment reference with a common MSbar coupling input.
Version 0.11 carries those channels and conventions into the full coupled
fiducial hard convolution in this one-family unpolarised-vector reference.
It restores the entire O(alpha) overlap and validates its second expansion,
causal projection, signed FF quadrature and scale/grid/sliver refinements.
See [the fiducial matching contract](docs/coupled_fiducial_matching.md).
Version 0.12 adds the matched pair-mass spectrum with a shared cumulative
measurement, conservative resolved-cell integration and correlated
fragmentation/starting-scale scans. Bin-specific refinements expose numerical
effects that cancel in the total rate; see the
[mass-spectrum benchmark](benchmarks/mass_spectrum_reference.md).
Version 0.13 adds both-muon pT and signed eta spectra with the same pair
selection, full hard matching and independent bin/geometry checks; see the
[single-muon benchmark](benchmarks/lepton_spectrum_reference.md).
Version 0.14 stabilises the endpoint panel action by evaluating endpoint strips
on the FF lattice and aligning shared panels with pT/eta contacts. The geometry
in C2 is now panel independent; see the [numerical contract](docs/endpoint_quadrature.md)
and [benchmark](benchmarks/endpoint_quadrature_reference.md).
The fixed-grid four-to-eight-panel maximum pT shift is now 3.25e-8 versus
1.16e-4 previously; the last grid doubling still changes a bin by 1.11e-5.
The next work is the remaining grid/interior numerical budget, a fuller
physical uncertainty budget, spin/flavour/threshold content and a demonstrated
production-host interface.
Develop and validate this small calculation before modifying a full NNLO program.

### Delivery stages and accuracy claims

| Stage | Concrete deliverable | Claim and limitation |
| --- | --- | --- |
| A: reference | Fixed-coupling, photonic/non-singlet FF evolution; Born Z decay; distribution-aware convolution | LL and fixed-order checks in a restricted QED model; no full-QED NLL claim |
| B: matched decay | Time-like NLL FFs and subtracted O(α) Z decay, initially in a one-lepton-family QED model | O(α) decay matched to NLL collinear resummation in the stated theory and pole approximation |
| C: phenomenological FFs | Explicit active-flavour, coupling, threshold and extra-pair prescriptions; required channels validated | NLL only for the documented physical content; unresolved pair/hadronic effects remain explicit limitations |
| D0: NLO W/Z host | Working W±/Z FF-only laboratory pilot and common-production PHOTOS comparisons; finite hard matches next | Matched NLL requires W/Z hard decay and physical-content gates; no full NLO EW claim |
| D: NNLO host | NNLO QCD production coupled to the validated decay/FF layer | NNLO QCD plus mixed logarithms demonstrated by expansion; NLL for the production–decay approximation only where its hard coefficients are supplied |
| E: extensions | Full off-shell NLO EW matching and, if accessible, an exact mixed correction | Exact fixed-order coverage listed by perturbative order; multiplying correction factors does not supply full NNLO QCD × NLO EW |

Stages A–B provide the existing restricted-model references. D0 is the next hadronic milestone; its FF-only host pilot can proceed before C, while a physical full-QED NLL claim still requires C and the W-specific matching checks. D follows validated NLO integration. Alternative schemes, full off-shell EW, exact mixed corrections and experimental dressing studies have separate gates.

## Physics contract to fix before coding kernels

### Conventions and observable definition

Use `D[parent -> observed](z, mu_frag)` throughout. Record the mapping to each paper's indices: PDFs and FFs use different evolution orientations. Use `a = alpha/(2*pi)` in algebra and avoid reusing this symbol for a parent index.

Write the convolution as a functional of a measurement `F`, avoiding unspecified Jacobians in a differential formula:

```latex
\sigma[F] = \sum_{i,j}\int d\Phi\,d\hat\sigma_{ij}(\Phi;\mu_{\rm frag})
\int_0^1 dz_1 dz_2\,D_{i\to\mu^+}(z_1;\mu_{\rm frag})
D_{j\to\mu^-}(z_2;\mu_{\rm frag})
F(z_1 k_i,z_2 k_j,X).
```

The host phase space and weights describe parent momenta `k`; `F` applies cuts and fills bins using observed momenta `p = z*k`. Any density in observed variables needs its own change-of-variable Jacobian; do not insert an extra `1/z` into this event-functional convention. At leading power, the two FFs describe two distinct hard directions. They do not specify every multi-lepton selection or two observed leptons arising from one parent.

Baseline: √s = 13 TeV, opposite-charge bare muons, each pT > 25 GeV and |η| < 2.4, and 60 < mμμ < 120 GeV. Start without photon dressing, isolation, photon vetoes or extra-lepton vetoes. Define bin edges and leading/subleading-muon ordering in configuration. Use an inclusive decay benchmark without hadronic PDFs for debugging.

The host must expose parent phase space **before observed-muon cuts**. Radiative loss moves events between bins, so already-cut histograms or accepted events with unknown generation cuts are insufficient. Hard real photons remain explicit in the matched calculation; collinear FF radiation is inclusive unless an observable extension supplies more information.

### Photonic, singlet and pair content

The first test may omit fermion loops and real pairs, using fixed α in a photonic/non-singlet model. This isolates distribution algebra, double convolution and matching without claiming full physical accuracy.

For the first NLL coupled calculation use a closed one-family basis `{μ⁻, μ⁺, γ}` and consistent running coupling. Store the full parent/daughter matrix, including photon survival and both charges, even if only muon columns enter the first observable. Extend to the physical lepton content afterwards.

Before promoting singlet evolution to a fiducial prediction, specify how additional μ⁺μ⁻ pairs are counted. An inclusive identified-lepton density counts multiplicities; it is not automatically an exclusive event with exactly two muons or the pair closest to MZ. Establish by order counting which hard photon/other-parent channels enter the claimed towers. If selection correlates particles from one collinear sector, assess a multi-particle fragmentation or resolved-pair treatment. A tagged-primary/non-singlet result is an allowed interim deliverable, labelled accordingly.

Choose a documented prescription for e, μ and τ thresholds and light-quark vacuum polarisation. Perturbative light-quark thresholds around ΛQCD do not replace the physical hadronic contribution to α. Distinguish species in the evolution matrix from those in the beta function and justify consistency. Including quark/gluon fragmentation requires a broader QCD+QED system and boundary conditions; adding quarks to running α alone does not implement those channels. Joint-factorisation papers provide context, not proof that the low-scale boundary is purely perturbative.

### Coupling, scales and factorisation scheme

Use separate names for QCD renormalisation `mu_R`, proton factorisation `mu_pdf`, final-state fragmentation `mu_frag` and FF starting scale `mu_0`. Initially set central `mu_frag = MZ`. A later dynamic scale should depend on parent hard kinematics rather than post-fragmentation mμμ unless explicitly derived.

Use Gμ Born normalisation and α(0) for the reference O(α) photonic decay. For running-coupling NLL evolution, specify the renormalised coupling, boundary, beta coefficients, threshold matching and finite conversion to the hard calculation. Combining an α(0) hard coefficient with running-α FFs requires checking that conversion at the claimed accuracy.

Start in **MS̄**. The Δ-scheme of [Frixione 2021](https://arxiv.org/abs/2105.06688) is a PDF construction. A time-like extension is a derivation task: transform FF boundaries, NLO kernels and hard coefficients together, then show agreement of physical observables through the retained order. Endpoint behaviour of an unphysical FF alone does not select a scheme. Defer the extension until MS̄ matching passes.

## Matching and perturbative bookkeeping

### O(α) decay matching

Let `H0` be the Born decay, `D = I + a D1 + ...`, and `H1_sub` the coefficient of the massless, collinearly subtracted O(α) hard decay. Then

```latex
R_{\rm dec} = (H_0 + a H_{1,\rm sub})\otimes D\otimes D,
\qquad
[R_{\rm dec}]_{a^1}
= a H_{1,\rm sub}
+ a H_0\otimes(D_1\otimes I+I\otimes D_1).
```

Use the same measurement and mappings in both terms. Their sum must reproduce the leading-power massive O(α) decay, including finite terms, not just mass logarithms. Massive additive matching is also possible, with overlap defined by expansion of this same resummed approximation. Do not attach an additional EW FSR shower without its own overlap prescription.

### QCD orders and optional exact mixed matching

Represent hard coefficients by `H[k,n]`, relative orders αs^k α^n, keeping the Born EW prefactor separate. The first NLO QCD W host supplies `H[0,0]` and `H[1,0]`; a later NNLO host adds `H[2,0]`. FF convolution generates universal mixed collinear terms, but a full NLL cross-section claim at each QCD order also needs the appropriate O(α) hard information: `H[k,1]` feeds higher-order NLL towers. A pole-approximation production–decay construction can supply an approximation to these coefficients; a plain QCD event file cannot. The W pilot has one observed charged leg rather than the two-leg neutral-current action above.

For an off-shell fixed-order anchor use an explicit additive construction. With `R` the chosen resummed prediction already containing pure NNLO QCD, for example:

```latex
\sigma_{\rm matched} = R
+ \big(\Delta\sigma_{\rm NLO\,EW}^{\rm exact}-[R]_{(0,1)}\big)
+ \big(\Delta\sigma_{\rm mixed}^{\rm exact}-[R]_{(1,1)}\big).
```

The last bracket is optional. Each Δσ is an order-specific correction, not a total cross section. Verify coefficients through `(2,0)`, `(0,1)` and, if included, `(1,1)`. This fixes those orders; it does not upgrade every higher-order tower to exact off-shell NLL.

Subtract the **entire** `(1,1)` coefficient of `R`: both muon legs, finite O(α) FF boundaries, hard-decay terms, coupling conversions and included channels. Subtracting only `Pμμ * log(mu_frag²/mμ²)` leaves overlap. The exact-minus-approximate residual is a **mixed remainder**, containing factorizable finite effects as well as possible non-factorisable contributions; do not call it purely non-factorisable. Align PDFs, EW inputs, cuts and normalisations.

Make order extraction a first-class operation. Coupling scans are diagnostics; analytic/coefficient-level expansion is the acceptance reference.

## WP0 — source audit, backend spike and host discovery

Deliver a conventions note and dependency/capability table with pinned revisions, plus an executable host probe when access exists.

1. Transcribe O(α) FF boundaries from [Frixione 2019, section 6](https://arxiv.org/abs/1909.03886), recording equation numbers, plus prescriptions and index mapping. Check flavour/mass assumptions before replacing me by mμ. Keep formula provenance beside the implementation.
2. Locate explicit time-like NLO QED kernels or derive them by diagram-level abelianisation. A global `CF -> 1, CA -> 0, TR*nf -> sum Qf²` substitution does not document all off-diagonal, pure-singlet and charge factors. Produce a channel-by-channel source table and independent moment checks.
3. Compare eMELA reuse with a small standalone implementation. [eMELA](https://github.com/gstagnit/eMELA) documents space-like PDFs, Fortran evolution and a C++ endpoint wrapper. Time-like kernels and endpoint solutions must be implemented/derived; they are not obtained by toggling a flag. Audit build, licence and extension points before choosing a fork.
4. Treat [APFEL++](https://github.com/vbertone/apfelxx) as another candidate. Its [documented initialisers](https://vbertone.github.io/apfelxx/html/namespaceapfel.html) include time-like QCD; that does not establish a ready-made time-like NLO QED FF backend. Verify the required API and kernel content in a pinned revision.
5. First probe local MCFM W± at NLO QCD for an analysis hook with signed real/subtraction contributions, then consider MATRIX, NNLOjet or a public alternative for NNLO. Require a small run reproducing native inclusive/fiducial results. Record access, licence, modes, spin information, generation cuts, slicing/subtraction controls and photon content. Do not assume every code exports reusable events.

Prefer C++17/CMake for the core and host interface, with Python for analytic references, scans and plots. Keep any Fortran dependency behind the evolution adapter. Choose one production backend after the spike; retain a small independent fixed-order/moment reference rather than implementing two full backends initially.

## WP1 — distribution algebra, boundaries and evolution

### WP1a: distributions and reference calculation

An NLO FF contains δ(1−z), plus distributions and possibly signed pieces. It is not generally a positive probability density. Support distribution action and moments, not only `D(z, mu)`.

Implement regular, delta and plus pieces with

```latex
\int_0^1 dz\,[g(z)]_+ f(z)
=\int_0^1 dz\,g(z)\,[f(z)-f(1)].
```

Include the endpoint contribution for finite lower bounds. Resolve cut/bin discontinuities using breakpoints or adaptive integration. Use analytic endpoint subtraction and common quadrature points for correlated terms. For later Monte Carlo sampling, use a positive proposal with signed weights and explicit endpoint contributions. Implement the two-leg subtraction consistently, including its double-endpoint term.

First use an isotropic Born decay and analytic cut examples. Required leading-power kinematic checks: `pT' = z*pT`, η unchanged for each muon, and `mμμ'^2 = z1*z2*m_parent-pair²`. The isotropic source is a numerical fixture; physical production later requires spin correlations.

### WP1b: NLL evolution

Reference progress: the one-family matrix, running real moments and O(α²)
coefficient extraction are implemented in version 0.8; see the
[solver contract](docs/moment_evolution.md). Version 0.9 adds the
[x-space reference](docs/xspace_evolution.md), endpoint projection, moment
recovery and independent causal checks. Broader measurement and full
hard-decay matching below remain required.

NLL requires O(α) boundaries, LO and NLO **time-like** kernels, and running/matching terms appropriate to the chosen counting, including the two-loop QED beta function in the standard running-coupling NLL solution. Retain fixed coupling as an independent diagnostic. [Bertone et al.](https://arxiv.org/abs/1911.12040) provide space-like NLL methodology, not replacement time-like kernels.

Choose Mellin-space evolution with controlled inversion or an x-space method after WP0. Derive the time-like large-z form or an equivalent stable endpoint treatment. A grid excluding z=1 without its unresolved contribution is insufficient. Separate numerical endpoint cutoffs from physical parameters and demonstrate convergence/cancellation.

Expose at least:

```text
initial_condition(parent, observed, mu_0, config)
moment(parent, observed, N, mu_frag, config)
apply(parent, observed, measurement, mu_frag, config)
expand(alpha_order, parent, observed, mu_frag, config)
```

Configuration contains theory content, scheme, coupling prescription, masses, thresholds and numerical tolerances. Persist it with grids; reject incompatible grid/hard-coefficient combinations.

### WP1 validation

- Reproduce boundaries and O(α) scale dependence. Expand to O(α²) too: verify LL `L²` and NLL `L` coefficients, including boundary × evolution and beta terms. An O(α) check cannot test an NLO kernel.
- In a scalar non-singlet fixed-coupling diagnostic, verify `D = I + a(d1 + L P0) + a²[L P1 + L P0⊗d1 + L² P0⊗P0/2]` through the displayed logarithms, with `L = log(mu_frag²/mu_0²)`. The unknown two-loop boundary constant is beyond NLL. Include the additional running terms for running-coupling checks.
- Check evolution composition, equal-scale identity and threshold matching in moments. Compare x-space actions/inversion with an independent reference at representative scales and near the endpoint.
- Check charge/valence conservation at N=1 and momentum conservation at N=2: for **each parent**, sum `∫ dz z D[parent -> h] = 1` over the complete implemented daughter set. A single observed-muon FF does not satisfy that momentum sum rule by itself.
- Target moment residuals below 10⁻⁶ in the closed reference theory, supported by integration/step refinement. Do not use an incomplete daughter sum as a physical conservation check.
- Compare LL logarithmic coefficients with [Dittmaier–Huber](https://arxiv.org/abs/0911.2329), aligning coupling, soft-exponentiation and finite-term conventions. Complete formulas need not agree when they include different beyond-LL terms.

## WP2 — subtracted and massive O(α) decay benchmark

Deliver massless subtracted Z → μμ(γ), an independent massive O(α) decay reference, and an inclusive/fiducial validation report.

1. Begin with photonic decay in the pole approximation. Preserve the Z spin-density matrix for the eventual production–decay interface. Fix width/branching-ratio normalisation and expand normalised decay factors consistently; avoid applying the same EW decay correction twice.
2. Choose one identified-lepton subtraction method, preferably compatible with the eventual host. FKS/identified-particle dipoles are candidates; [fragmentation antennae](https://arxiv.org/abs/2201.06982) provide another route. Document local mappings, integrated poles, finite remnants and FF counterterms. Fix signs and MS̄ normalisation from the complete calculation rather than a schematic `P/epsilon` term.
3. Verify soft/collinear cancellation and subtraction-parameter independence. Use the same measurement contract for real and mapped counterevents. The first analytic photonic decay can supply its virtual term directly without a full EW virtual provider.
4. Implement the massive reference with identical inputs, or establish code access permitting a genuine lepton-mass scan. Demonstrate the required FEWZ/MATRIX mode instead of inferring it from published numbers. Compare identical pole approximations first; an off-shell-versus-pole comparison is not a mass-power test.
5. Scan `{mμ, 0.25, 0.5, 1, 2} GeV`, refining the range as needed. Fit power corrections allowing `m² log^k(m²/Q²)` terms. Cuts/endpoints can amplify power effects; do not require a pure `C*m²` fit in every narrow bin. Retain numerical covariance and reference-code systematics.
6. Show O(α) cancellation of `mu_frag` dependence between hard and FF terms. Vary `mu_0` together with its logarithmic boundary and coupling matching, checking starting-scale cancellation through the retained order.

Target relative agreement at physical mμ: 10⁻⁴ for the fiducial rate and selected populated bins, once combined numerical errors are smaller. Use absolute errors or a stated reference normalisation in nearly empty bins. Report achieved precision if the independent reference cannot resolve this target.

Full off-shell NLO EW is a separate extension using OpenLoops/Recola or an existing provider: include weak terms, ISR, IFI and required photon-induced channels. Match initial-state subtraction to the **specified proton PDF scheme**; a LUX-derived photon density alone does not specify every counterterm convention. Keep initial-state and final-state scales/schemes separate.

## WP3 — NLO W production adapter, then NNLO extension

### Host contract

Start with the [NLO W milestone](docs/w_nlo_milestone.md): local MCFM 10.3.1 is the first candidate to exercise. Its [W process documentation](https://mcfm.fnal.gov/processes/process1.html) specifies NLO support, including dipole subtraction. A working adapter and identity test remain required. For the later NNLO extension, MATRIX and NNLOjet are candidates; [the MATRIX mixed calculation](https://arxiv.org/abs/2412.16095) establishes a physics reference, not local access or an export API.

The internal record or callback must retain:

```text
parent four-momenta/flavours; explicit hard photons and QCD final state
signed weights, units, perturbative-order labels, scale/PDF metadata
correlation groups for real/subtraction/projection contributions
spin/decay information when attaching an independently generated decay
generation cuts, subtraction/slicing settings, configuration identifier
```

Choose the external format after inspecting the host. Do not require unweighted events or HepMC: NNLO calculations often expose signed, correlated contributions. Prefer an in-process analysis hook when export loses cancellation information or costs too much.

For post-processing, demonstrate sufficient differential information and all required contributions. An LHE file or Born reweighting is not automatically NNLO accurate for arbitrary post-fragmentation observables. Pure-QCD lepton momenta/weights support a restricted FF-convolution result; attaching a spin-correlated hard decay needs additional information or host integration.

### Integration sequence

1. With `D = I`, reproduce native inclusive/fiducial results, order by order where possible. Verify units, channels, negative weights and scale variations.
2. Add expanded O(α) FF action before cuts. Keep counterevents and variations correlated; use shared z nodes/random variables where mappings permit. Measure cancellation quality and variance before production.
3. Add resummed action through deterministic integration/tabulation first. Profile endpoint integration, interpolation and sample counts; added cost must be measured rather than assumed negligible.
4. Attach the validated pole decay with production spin correlations and QCD recoil. Verify its expansion/normalisation and list the `H[k,n]` coefficients actually supplied.
5. Check host-regulator convergence after fragmentation: qT slicing/extrapolation if applicable, subtraction parameters, generation cuts and FF tolerances. Symmetric pT cuts and Jacobian regions need particular attention. Use mildly asymmetric pT cuts as a diagnostic while keeping the nominal definition explicit.
6. If an exact mixed correction and suitable mass scans are accessible, compare the leading mass-log coefficient with the `(1,1)` expansion, then apply complete-overlap matching. Diagnose conventions, channels, mappings and regulators before suggesting factorisation failure.

Baseline cuts on observed collinear muon momenta do not themselves invalidate leading-power factorisation. Isolation, resolved photons, very small pair pT and other recoil-sensitive measurements require their own measurement/factorisation analysis; exclude them from the initial accuracy claim.

## Software layout and reproducibility

Create only components needed for the current stage. Proposed layout:

```text
CMakeLists.txt
include/leptonff/              # config, distributions, evolution, expansion
src/                          # core and chosen backend adapter
apps/                         # moments, decay benchmark, host-analysis drivers
python/                       # analytic references, scans, plots
configs/                      # versioned reference/fiducial setups
tests/                        # analytic and integration checks
benchmarks/                   # compact reference inputs/results and provenance
docs/                         # conventions, kernel sources, accuracy, host contract
```

Core distribution/evolution checks must run without NNLO codes, LHAPDF or EW providers. Keep analytic references independent of production code where practical.

Record configuration, source/dependency revisions, compiler, backend, seeds/quadrature, cuts/bins, units, achieved errors and host regulator settings with every result. Preserve covariance between central/varied predictions. Large grids/event files are generated artifacts with manifests, not unexplained fixtures.

Document commands for building, fast checks, the decay benchmark and one reproducible plot/table. Keep long host integrations separate from fast regression checks.

## WP4 — uncertainty, comparisons and physics output

Vary the **matched hard and FF pieces together** with `mu_frag`; varying only the FF leaves uncancelled dependence at a retained order. Start with `mu_frag/MZ = {0.5, 1, 2}`, independently of standard QCD `mu_R`/`mu_pdf` variations. For `mu_0/mμ = {0.5, 1, 2}`, recompute boundaries and threshold/coupling matching consistently; use the lower point only where the prescription defines it, otherwise document an admissible range.

Report numerical errors, QCD-scale variation, fragmentation-scale variation and starting-scale sensitivity separately. Assess pole/off-shell, pair/hadronic and missing wide-angle-photon effects where relevant. A small fragmentation band does not bound all omitted effects. Wide-angle soft radiation is not generally suppressed by mμ²/Q² in a fiducial cross section.

| Comparison | Purpose | Acceptance or interpretation |
| --- | --- | --- |
| Moments and O(α²) expansion | FF correctness and NLL content | Analytic coefficients, conservation identities, numerical convergence |
| Massive vs expanded matched decay | Finite terms, scheme, mapping | Leading-power agreement with resolved numerical/power effects |
| Native host vs identity convolution | Adapter correctness | Agreement within correlated integration errors |
| Mixed mass-log coefficient, if accessible | QCD × collinear overlap | Coefficient agreement with aligned conventions; optional for a labelled restricted result |
| LL vs NLL with identical hard matching | Added collinear accuracy | Resolve differences where possible; do not presume size/sign |
| Photos and Sherpa YFS with documented settings | Resummation/matching comparison | Same production, cuts, pair content and EW inputs; avoid duplicate FSR |

Photos and YFS implementations have process-dependent matrix-element matching and soft-resummation content. The motivation is a transparent time-like collinear expansion and matching study, not a blanket claim that these tools lack perturbative accuracy. [Sherpa's decay study](https://link.springer.com/article/10.1140/epjc/s10052-019-6614-x) explicitly includes higher fixed-order decay corrections.

Deliver rate and mμμ, pTμ, ημ ratios for fixed order, matched LL and matched NLL, with an accuracy table and achieved error budgets. Mark unstable endpoint bins. The log `log(MZ²/mμ²) ≈ 13.5` motivates resummation; cut-dependent shifts and scale bands must be calculated, not installed as expected results.

Bare-to-dressed ratios need explicit photons, dressing radii and a validated measurement extension. Longitudinal FF rescaling alone cannot predict arbitrary cone dependence. Fixed nonzero-radius dressing is collinear safe to all orders with respect to lepton–photon collinear radiation; extra-pair flavour definitions may require further care. Defer this study.

## Gates, dependencies and schedule

| Gate | Evidence required | If unresolved |
| --- | --- | --- |
| G0 conventions/backend | Kernel provenance, indices/coupling definitions, closed reference theory, working backend spike | Continue analytic/distribution work without an unverified NLL claim |
| G1 FFs | O(α²) LL/NLL coefficients, moments, endpoint/step convergence | Fix kernels, boundaries or numerics before matched production |
| G2 matched decay | Massive/expanded agreement, finite terms, regulator stability, O(α) scale cancellation | Debug in MS̄ and the simplest photonic theory before adding complexity |
| G2b physical content | Flavour, hadronic-coupling and pair-selection prescription; channel/order counting | Release a restricted-model/tagged-primary result with restricted accuracy |
| GH0 NLO W host | Identity and expanded one-leg convolution reproduce/check native NLO results with correlated contributions, cuts and uncertainties | Retain FF-only pilot label until W hard matching and physical content are validated |
| GH host | Extend the validated NLO adapter to a working NNLO hook/output with correlated convolution/regulator checks | Retain validated lower-order deliverable before claiming NNLO DY |
| G3 exact mixed | Exact correction accessible; entire `(1,1)` overlap removed | Omit exact mixed matching and document approximate terms retained |
| G4 release | Reproducible commands, accuracy table, resolved numerical errors, omitted effects documented | Narrow claims or gather statistics; never require a preferred effect size |

Next dependency chain: NLO W host/identity probe → one-leg FF-only pilot → W hard-decay matching plus physical-content prescription → aligned W comparisons → NNLO extension. The W hard-decay work can proceed alongside the pilot. Existing neutral-current grid refinement continues independently. Alternative schemes/exact mixed matching do not block the first validated decay release.

Indicative scientific effort, to be re-estimated after backend/host probes:

- WP0 and WP1a: 2–4 weeks for sources, capability audit, distributions and toy decay.
- WP1b and reference-model NLL validation: 1–3 months, depending on kernels/endpoints.
- WP2 and physical-content prescription: 2–4 months, with independent massive validation as a dependency.
- WP3–WP4: 2–6 months after a usable host interface, including statistics/regulator scans.
- Full off-shell, further physical singlet/pair extensions and exact mixed matching: estimate separately after access/theory requirements are established.

These are research effort ranges, not a guaranteed delivery date. Kernel derivation, endpoint control, the independent massive reference, physical pair content and host access carry material risk; none is dismissed as pure bookkeeping.

## First implementation batch

Steps 1–4 are available as a fixed-order reference, followed by standalone
fixed-coupling non-singlet LL evolution with grid/cutoff controls. See
[README.md](README.md), [implementation status](docs/status.md), and the
[reference-backend decision](docs/backend_decision.md). Full fiducial NLL
hard/pair matching, physical spin-correlated decay and host probes remain
outstanding.
The evolved O(α) diagonal boundary and complete FF-only O(α) overlap matching
are also implemented; [details](docs/boundary_matching.md) and the
[NLO kernel audit](docs/nlo_kernel_audit.md) identify the remaining NLL ingredients.
The [mass-only vector hard benchmark](docs/vector_mass_matching.md) now supplies
finite hard+FF scale cancellation, a massive real/virtual reference and a genuine
mass scan. The [fiducial vector extension](docs/fiducial_vector_matching.md)
adds symmetric angular cuts, an exact massive reference, and additive NLO+LL
matching with full O(alpha) overlap checks. The [hard-dressing extension](docs/hard_dressing.md)
now convolves H1 with both LL-evolved muons, with independent signed/Poisson
and angular checks. This is still restricted-model G2 evidence: the observable
NLL evolution/matching and physical spin/current and boosted-parent extensions
remain required.
The [NLO lepton-block increment](docs/nlo_kernel_audit.md) now supplies explicit
same/opposite-charge kernels and cancels the singular singlet analytically in
the conserving valence channel. The source-separated photonic first moment
does not vanish; do not force its normalisation or insert it into the scalar
LL solver. Version 0.8 completes the [photon channels](docs/nlo_photon_kernels.md),
per-parent momentum checks and [coupled real-moment evolution](docs/moment_evolution.md)
with two-loop MS̄ running and explicit O(α²) extraction. Version 0.9 supplies
the endpoint-controlled coupled x-space action, moment recovery and a Born
two-direction FF fiducial match with independent causal/grid/domain checks.
Version 0.10 reconstructs the hard photon-primary coefficients and subtraction,
validates their fixed-fraction fiducial action, and combines all hard-parent
channels with coupled FFs for inclusive pair-mass moments. The hard coupling
uses the same MSbar input as the FFs. Pair counting is inclusive across distinct
hard parents, with no extra-pair veto. Next, carry this contract into the full
fiducial x-space hard convolution and validate the combined numerical
refinements. Version 0.11 now delivers that restricted-model fiducial rate:
all hard-parent channels, common MSbar coupling, complete first overlap and
explicit second-order extraction pass independent causal/geometry checks.
Its numerical refinements are recorded in
[the benchmark](benchmarks/coupled_decay_reference.md). The physical Z spin,
multi-flavour and NNLO-host milestones remain open.
Version 0.12 extends this reference to matched pair-mass bins with common
cumulative quadrature, conservative resolved-cell integration, full component
and expansion closure, and per-bin numerical/fragmentation/starting-scale
studies. Version 0.13 adds the single-muon pT and signed eta spectra with the
same pair selection and two entries per inclusive pair. Selected one-lepton
cumulatives share all source/projection data and carry the complete hard/FF
match. Literal signed lattice, allowed-azimuth and selected massive checks
accompany bin-specific numerical and correlated scale scans; see the
[contract](docs/lepton_spectra.md) and [benchmark](benchmarks/lepton_spectrum_reference.md).
Version 0.14 adds direct endpoint geometry and measurement-aligned common
panels. It removes the leading endpoint interpolation error while retaining
signed linear projections and the complete numerical expansion. For every
spectrum, bin numerical shifts must not be replaced by the smaller total-rate
shift. Physical flavours/spin/boosts, a complete uncertainty budget,
external comparison and a validated production-host adapter remain open.

1. Create C++/CMake skeleton and conventions note: indices, `alpha/(2*pi)`, scale names.
2. Implement distribution action/moments with analytic fixtures for delta, plus kernels, finite lower bounds and two-leg cuts.
3. Transcribe MS̄ one-family O(α) boundaries with source equations; verify charge/momentum moments.
4. Build Born-decay driver demonstrating identity convolution and pT/mμμ rescaling, then add O(α) FF action.
5. Finish the backend spike and pin the chosen dependency before NLL implementation; explicitly list missing time-like NLO kernels.
6. Establish host and massive-reference probes as independent workstreams; neither capability is assumed until exercised.

This batch ends with a reviewable numerical reference and settled interfaces. It does not begin by modifying a private NNLO code or deriving a new scheme.

## Reference audit and reading map

Checked against primary pages/source descriptions on 17 October 2026. Bibliographic relevance is distinct from verified runtime capability; pin equations and software revisions during WP0.

| Reference | Use and qualification |
| --- | --- |
| [1909.03886 — Frixione](https://arxiv.org/abs/1909.03886) | O(α) FF boundaries, section 6; primary starting point |
| [hep-ph/0205172 — Arbuzov–Melnikov](https://arxiv.org/abs/hep-ph/0205172), [hep-ph/0206036v2 — Arbuzov](https://arxiv.org/abs/hep-ph/0206036v2) | Implemented NLO time-like charged-lepton block; pair-interference and valence checks; not a full photon/lepton matrix |
| [1606.02887v3 — de Florian et al.](https://arxiv.org/abs/1606.02887v3), [hep-ph/9612250v1 — Stratmann–Vogelsang](https://arxiv.org/abs/hep-ph/9612250v1) | Space-like QED plus full analytic continuation gives the implemented photon channels; independent time-like formula and moment checks |
| [1911.12040 — Bertone et al.](https://arxiv.org/abs/1911.12040) | Space-like NLL numerical/analytic methodology, not time-like kernels |
| [2105.06688 — Frixione](https://arxiv.org/abs/2105.06688) | PDF Δ-scheme; time-like counterpart needs consistent derivation |
| [eMELA](https://github.com/gstagnit/eMELA), [APFEL++](https://github.com/vbertone/apfelxx) | Candidate infrastructure; required FF capabilities need verification |
| [0911.2329 — Dittmaier–Huber](https://arxiv.org/abs/0911.2329), [1511.08016 — Dittmaier–Huss–Schwinn](https://arxiv.org/abs/1511.08016) | FSR structure-function and resonance-region mixed benchmarks |
| [2201.06982 — Gehrmann–Schürmann](https://arxiv.org/abs/2201.06982) | Identified-photon fragmentation subtraction; QED adaptation is work |
| [1208.5967 — FEWZ](https://arxiv.org/abs/1208.5967) | EW/QCD combination reference; verify required mass scan/export separately |
| [2412.16095 — neutral-current mixed calculation](https://arxiv.org/abs/2412.16095) | Massive-lepton mixed reference; access/input flexibility remain dependencies |
| [2201.01754 — Armadillo et al.](https://arxiv.org/abs/2201.01754) | Neutral-current two-loop mixed virtuals; replaces 2405.00612, which is charged-current DY |
| [2505.23487](https://arxiv.org/abs/2505.23487), [2607.07664](https://arxiv.org/abs/2607.07664) | Joint QED+QCD factorisation context. The latter concerns single-inclusive hadron/jet production in lepton–hadron scattering and uses LO evolution in its numerical analysis; not an established independent NLL time-like muon-FF benchmark |
| [2508.16964](https://arxiv.org/abs/2508.16964) | NNLO electron PDFs from SCET; PDF boundaries are not available NNLO FF boundaries |
| [Sherpa decay study](https://link.springer.com/article/10.1140/epjc/s10052-019-6614-x) | Soft-resummed comparison with specified fixed-order options |
| [2310.11574 — ATLAS Z → ℓℓγ](https://arxiv.org/abs/2310.11574) | Experimental FSR context; not itself a supplied unfolding-factor table |

Novelty and publishability remain literature-review questions. Do not retain unverified claims that no comparable FF code exists, that this is the first fixed-accuracy FSR statement, or that particular generator differences/scale bands are already known.

## Version 0.18 coupled physical-current increment

The W± and charge-symmetric Z/gamma* spin operators now act on the complete
one-family coupled evolution: all hard-parent channels, inclusive pair
counting, the common MSbar coupling at original M, continuum O(alpha) overlap
and explicit C2. See [the contract](docs/coupled_lab_matching.md) and
[benchmark](benchmarks/coupled_lab_reference.md). This closes the split between
the earlier physical-current photonic pilot and the restricted one-family
NLL backend. It does not close physical flavour thresholds/scheme inputs,
experimental pair selection, hard differential spectra or a numerical and
statistical phenomenology budget. NLO QCD host validation remains the first
production milestone; NNLO QCD follows separately.

The completed v0.18 audit passes 52 core tests, 486 independent polarised
RECOLA points and 30 host controls (including 8192/16384-bin points and
mu_frag=M/2,2M). NLO QCD W± and Z matching is operational for the stated
one-family inclusive observables. The next implementation priorities are:

1. Reduce the hard-estimator variance and control the measurement/grid limit;
   current higher-order shifts are comparable to numerical/statistical changes.
2. Add hard-event differential pT/eta/mT and pair-mass measurements, retaining
   common random streams and bin covariance rather than labelling original
   Born directions as radiative spectra.
3. Extend active flavours/thresholds and fix the physical coupling input and
   pair-selection prescription before a quantitative data comparison.
4. Build the joint QCD/fragmentation-scale, numerical and statistical budget;
   retain the fixed-alpha photonic branch for the aligned PHOTOS benchmark.

A narrow fragmentation-scale band is not a substitute for unresolved grid
or rare-event errors. The one-family coupled PHOTOS differences are diagnostic
because its pair content and running-coupling convention differ.

## Version 0.19 numerical measurement increment

The first numerical-precision item now has conservative signed resolved
cells, exact endpoint atoms, analytic two-leg mass/pT cut integration and
consistent C1/C2 restoration. See [contract](docs/conservative_fraction_cells.md)
and [benchmark](benchmarks/fraction_cell_reference.md). Default atomic mode
reproduces the previous rates. Two-seed W/Z variance studies and grids through
16384 bins accompany fixed-decay responses through 32768 bins. The cell
measure is a different finite-grid reconstruction, with the same intended
continuum limit; observed gains do not certify a global error bound.

Remaining numerical work includes rare-tail/angular variance and a full
measurement/grid limit. The separate functional milestones remain hard-event
differential spectra, physical flavour/threshold/input matching, and a joint
phenomenology uncertainty budget. No small scale band substitutes for those
checks.

The v0.19 campaign is complete: 53 core tests and 36 host checks pass, with
no diagnostic flags over 27 runs. Observed higher-order variance gains over
two seeds are 1.09–2.52 (W+), 1.09–1.30 (W-) and 1.18–1.20 (Z). Fixed-decay
cell responses also have smoother 4096–32768-bin convergence. This completes
the conservative-cell increment, not the full numerical/statistical budget;
rare-event and endpoint/continuum-rate precision remain open.


## Version 0.20 differential benchmark completed

The NLO QCD W±/Z matching milestone now includes physical observed-momentum
spectra, full grouped-event cross-bin/channel covariance and normalised shapes.
W has muon pT, absolute eta and mT. Z has charge-symmetric pair mass, minimum
muon pT and maximum absolute eta; these are not single-muon or arbitrary
experimental leading-pair measurements. Both the coupled one-family NLL and
the fixed-alpha photonic PHOTOS branches are implemented.

The 39-run coupled campaign passes 54 required checks; the 12-run aligned
photonic campaign passes 24. Nine high-chi-squared flags in the reduced-call
coupled campaign remain recorded. Larger central runs, seven-point QCD and
three-point fragmentation-scale studies, grid/cutoff/seed diagnostics and
standalone plots accompany the new event covariance. These close the
functional differential milestone, not precision certification. See the
[contract](docs/lab_spectra.md), [benchmark](benchmarks/lab_spectrum_reference.md)
and [operational summary](docs/operational_summary.md).

The implementation dependency chain is now:

1. Completed in v0.21: full charged-flavour evolution and hard matching.
   Separate produced-flavour multiplicities from closed-loop charge sums in
   every time-like P1 channel; validate flavour/charge and momentum sums before
   evolution. Changing only running alpha or applying one global nf factor is
   insufficient. Keep the one-family reference as a regression mode.
2. Leptonic thresholds and alpha-input matching completed in v0.21. Version
   0.22 adds sourced HVP matching and an explicit pointlike boundary, with
   optional LO QCD evolution. Mixed coefficients and fitted fragmentation
   concern the optional combined programme, not the pure-QED NLL target.
   Document which real pair channels
   are included and how they map to the requested observed lepton definition.
3. Establish the target differential precision: independent higher-statistics
   samples, grid/domain and hard/angular-cutoff limits, plus recoil treatment
   near Jacobian/low-qT regions. Retain event covariance and native-error
   reconstruction throughout.
4. Complete PDF/alpha_s uncertainties with the selected production PDF order;
   align energy, cuts, line shape and widths with a selected bare/post-FSR
   measurement. Add charge-summed single-muon Z spectra or parity-odd spin
   information when the requested observable needs it.
5. Extend the validated host to NNLO QCD separately. Full NLO EW and exact
   mixed corrections are not implied by final-state-QED fragmentation matching.

The immediate common-production analytic-versus-PHOTOS benchmark is usable
now through the photonic branch. A physical multi-flavour NLL data prediction
remains gated on the items above.


## Version 0.21 flavour and threshold implementation

The first dependency in the v0.20 chain is implemented: general charged-flavour
P0/P1 with distinct produced-colour and closed-loop factors, full-parent
moment/x-space evolution, finite time-like heavy/photon maps and matched
running alpha. The e/mu/tau preset, tau matching-scale variations, explicit
MSbar and leptonic on-shell inputs are integrated into the NLO QCD W±/Z host.
The hard alpha and C2 running logarithm use the same scheme and thresholds.
The previous one-family and photonic branches remain regression modes.

The second dependency is complete for leptonic inputs only. A physical
hadronic vacuum-polarisation input and QCD evolution of coloured parents still
need a separate implementation; custom coloured entries are pure-QED tests.
No perturbative light-quark pole masses are silently substituted for hadronic
vacuum polarisation. The remaining precision, extra-pair/data definitions,
PDF/alpha_s, recoil and NNLO host work in items 3–5 above still applies.
The validation passes 58 core tests and 90 required host checks across 36 runs;
twelve reduced-call W+ host convergence warnings remain recorded.
See [contract](docs/flavour_threshold_matching.md) and
[benchmark](benchmarks/flavour_threshold_reference.md) for verified scope.


## Version 0.22 hadronic-input and QCD increment

Implemented a sourced on-shell five-flavour HVP input and a fermionic MSbar
anchor, with a residual finite threshold map fitted after perturbative quark
running to avoid double counting. An independently compiled installed
Jegerlehner routine provides an alternative input. The expanded parent system
includes a gluon, LO time-like QCD quark/gluon kernels, variable-flavour one-loop
alpha_s running and QCD-dressed C1/C2 QED expansion coefficients. W±/Z host
controls retain full event covariance and the legacy regression modes.

The named `hadronic-pointlike` profile makes the zero low-scale light-quark
fragmentation boundary explicit. The sourced HVP input does not determine
that real-fragmentation boundary. Mixed alpha*alpha_s splitting/beta terms,
NLO QCD fragmentation matching and fitted nonperturbative input are still
open; do not mark the full physical combined-accuracy gate complete. The
precision/data-definition/PDF/recoil gates remain as above.
See [contract](docs/hadronic_qcd_matching.md) and
[benchmark](benchmarks/hadronic_qcd_reference.md).

The v0.22 validation closes with 60 passing core tests and 101 required host
checks over 39 control runs, retaining thirteen W+ convergence flags. Final
release W±/Z runs reproduce every rate, spectrum, raw iteration and event
cross-product bit for bit. These checks close the implementation increment
and preserve the explicit approximation and precision gates above.

## Scope correction and pure-QED differential comparison

The original target is NLO QCD production combined with NLL final-state QED,
not combined strong/QED fragmentation evolution. Earlier entries that made
mixed kernels and fitted nonperturbative quark fragmentation prerequisites
to this target are superseded by this scope correction. QED flavour
thresholds and electromagnetic input matching remain relevant. The v0.22
strong-evolution model remains an optional extension.

The [accuracy audit](docs/qed_accuracy_audit.md) states the source, boundary,
coupling, time-like threshold and hard-matching conventions. The active
campaign compares the aligned photonic branch to PHOTOS and separately shows
leptonic NLL and the explicit pointlike HVP prescription with QCD FF evolution
disabled. It increases the call budget tenfold, repeats independent seeds,
doubles grids, widens the evolution domain and tightens radiation cutoffs.
The deliverables are absolute and normalised W± pT/|eta| and Z pair-mass
figures, paired ratios/differences, covariance, bin tables, numerical
sensitivities and provenance. Measurement-specific cuts/pair definitions,
PDF uncertainties and precision recoil treatment remain data-use work.

The comparison is complete and documented in
[the benchmark](benchmarks/qed_comparison_reference.md): 144 main consistency
checks, 15 HVP-input checks and 48 bitwise adapter regressions pass. The domain
audit fixed a silently ignored photonic span setting. The W− variance study
led to ten independent seeds at a smaller technical QCD dipole parameter,
reducing the problematic ratio error from about 5.4% to 0.52%. Final figures,
normalised-shape covariance and bin tables are in
[the comparison report](results/qed_comparison_v023/comparison.md). Precision
limits are explicit: grid responses up to 0.6% in the coarser HVP Z control
and noisy geometric-cutoff responses do not support per-mille certification.
