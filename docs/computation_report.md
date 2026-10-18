# Computation ingredients and evidence report

Review date: 18 October 2026. This document describes the calculation that was
actually run, distinguishes it from optional library capabilities, and records
the limits of its validation. Large raw records stay in the local ignored
`results/` tree; small tracked benchmark snapshots and reproduction programs
remain in the repository. Nothing in this review constitutes an experimental
data fit or a complete electroweak prediction.

## Assessment

The project has a working generator-independent QED fragmentation and matched
decay library, working production adapters, and substantial analytic and
independent numerical checks. The normalized **first-order Z mass correction**
has a strong external massive-matrix-element check. The implementation of our
specified photonic resummation formula has a much more extensive internal audit
than before. Those are meaningful achievements.

The remaining distinction is between correctly implementing a specified
photonic NLL prescription and establishing that it is the appropriate NLL
description of a physical no-pair bare-muon observable. The KKMC comparison
resolves differences at tens to hundreds of ppm in normalized moments. Its
coupling/mass scans do not uniquely assign those differences to a perturbative
order. The completed McMule fixed-order benchmark now finds a second-order
remainder consistent with a constant in the mass logarithm, supporting the
implemented NLL coefficients for the tested inclusive vector moments. It does
not yet explain the complete resummed KKMC residual. Its results are recorded in
[the fixed-order benchmark](mcmule_fixed_order.md).

“Certainly correct” below means a stated identity or implementation has been
verified over its declared scope. It does not mean that a test suite can prove
an entire phenomenological approximation correct.

## The different calculations must be kept separate

| Calculation | Production | QED decay/radiation | Purpose |
|---|---|---|---|
| Current fine-bin Z comparison | Born QCD, common SANC production | Massive SANC NLO; PHOTOS; PYTHIA; leptonff O(alpha), matched LL, matched photonic NLL | Independent FSR comparison with a common baseline |
| Leptonff-only Z accuracy study | Same Born baseline | O(alpha), LL, NLL; fragmentation-scale and finite-completion variations | Resolve small increments and numerical effects |
| Fixed-Q radiator audit | Fixed parent Q; no PDF convolution | Normalized two-leg photonic radiator and KKMC FSR | Isolate radiation, coupling and mass dependence |
| New McMule check | Pure vector e+e- current at fixed Q | Fixed-order photonic muon-line corrections through alpha^2 | Separate fixed-order coefficients from resummation completions |
| DYTurbo interface studies | Native NLO/NNLO and NLO+NLL/NNLO+NNLL QCD | Leptonff/PHOTOS responses | Validate production coupling and reusable interfaces |
| Older MCFM studies | NLO QCD | Photonic or coupled lepton/flavour prescriptions | Signed production-group and differential-decay validation |

The recent precision Z plots are **not** DYTurbo NNLO+NNLL predictions. A PDF
set whose name contains `nnlo` does not upgrade Born production to NNLO.
“NLL” must always specify whether it refers to QCD recoil resummation or QED
mass-logarithm resummation. They are independent accuracy labels.

## Observable, inputs and normalization

The current hadronic benchmark uses pp collisions at 13 TeV,
`NNPDF31_nnlo_as_0118`, central member, and production factorization scale Q.
The fixed electromagnetic coupling is alpha=1/137.035999084 and the muon pole
mass is 0.1056583755 GeV. SANC electroweak Born parameters and widths are retained
in each native JSON record's `settings`; they are not silently replaced by
DYTurbo or KKMC defaults. Source versions and dependency hashes are retained in
the corresponding build/run manifests.

For the Z mass study, parent Q is restricted to 50–200 GeV. The measured bare
pair-mass window is 50–116 GeV, with 0.25 GeV visible bins and **no lepton pT or
eta cuts**. Photons are excluded from the muon momenta. There is no photon
recombination cone. Feed-in from parents between 116 and 200 GeV is included;
parents above 200 GeV are not. This is a declared validation domain, not an
established full physical mass-window prediction. Its extension must be tested
before comparison with data.

The W reference uses parent Q=60–200 GeV, muon pT>25 GeV and |eta|<2.4. Its
visible bins are 0.25 GeV in pT, 0.05 in |eta| and 0.5 GeV in transverse mass.
No missing-pT or transverse-mass cut is imposed. The neutrino transverse
momentum defines particle-level missing momentum. Fine W shape precision is
still an outstanding task; these configurations must not be mistaken for a
selected CMS measurement.

The decay correction uses a **supplied inclusive branching normalization**.
It redistributes events in lepton momentum; it does not predict an absolute
two-loop partial width. PHOTOS-like unitary normalization, normalized SANC
FSR, and leptonff shapes can therefore be compared. The SANC raw inclusive
O(alpha) rate correction is retained as a separate counterterm channel.
Restoring a rate correction requires a consistent physical branching/width
choice, rather than multiplying it into an already corrected branching ratio.

For a fixed parent define z=m_mumu^2/Q^2 and
R_N=integral dz z^(N-1) R(z). The normalization is R_1=1. Reported moment
differences in **ppm** mean an absolute difference times 10^6, unless explicitly
labelled relative. Thus -86 ppm in a moment near 0.964 is about -0.0089 percent
relative, not an 86-percent or 0.86-percent effect.

## Fragmentation and hard matching ingredients

### Factorization, boundary and leading kernel

The factorization convention is timelike MSbar, with parent-to-observed
indices, a=alpha/(2*pi), and natural boundary scale mu0=m_mu. The whole-plus
action is

    integral [g(z)]_+ F(z) dz = integral g(z) [F(z)-F(1)] dz.

The leading charged-lepton kernel is P0=[p(z)]_+, with
p(z)=(1+z^2)/(1-z). The diagonal one-loop boundary is

    d1(z,mu0) = [p(z)*(log(mu0^2/m_mu^2)-2*log(1-z)-1)]_+.

The finite boundary is retained. It is not sufficient to evolve a delta
boundary and call that the complete matched result. Spectral convolutions
have their usual dz/z measure. An event measurement applied directly to a
weighted hard event instead integrates D(z)*F(z*p); adding another 1/z to
that event action would double count a Jacobian.

The LL evolution resums powers of alpha*log(mu_frag^2/m_mu^2). In the photonic
benchmark alpha is fixed. Fermion-loop running and real pair channels are
excluded deliberately; a photonic NLL calculation is meaningful without them,
provided its sector definition and matching are physically consistent.

### Finite hard coefficient and two-leg matching

Leptonff contains a finite O(alpha) hard-decay response, including the
implemented spin-correlated W currents and parity-even, charge-symmetric
Z/vector response. The production/decay tensor basis and recoil convention are
part of the interface. For the inclusive invariant-mass reduction, the sum of
the one-leg boundary and projected hard coefficient at mu_frag=Q is

    C(z) = [p(z)*(log(z)-1)]_+.

The soft endpoint logarithms in the separate hard and boundary pieces cancel
in this sum. Comparisons use the **complete hard-matched observable**. Historical
`ff_order2` columns are FF-only auxiliaries and are not a reference for the
complete O(alpha^2) decay correction.

Across both measured muon legs, the default resummation retains one subleading
insertion in total: the O(alpha) boundary, the O(alpha) hard coefficient, or the
O(alpha^2*L) evolution insertion. It does not multiply two independently
truncated “NLL legs” and inadvertently retain their extra products. Hard real
radiation is dressed with LL evolution in this counting.

The O(alpha), O(alpha)+LL and O(alpha)+NLL curves are nested predictions, not
three independent corrections to add together. The NLL curve contains the LL
content plus its designated subleading increment. Coupled-flavour modes are
alternative physical prescriptions, not an extra correction to add on top of
the fixed-alpha photon-only curve.

### Photonic NLO kernel and conditioning

The source-separated timelike photonic distribution G follows
[Arbuzov–Melnikov, equation (15)](https://arxiv.org/html/hep-ph/0205172), with an
independent formula check against
[Arbuzov–Scherbakova, appendix A.2](https://arxiv.org/html/0706.2984).
For 0<z<1 its regular expression is

    G_regular = p(z)*(2*log(z)*log(1-z)-2*log(z)^2-2*Li2(1-z))
                +(1+z)*log(z)^2/2+2*z*log(z)-3*z+2.

Its source delta coefficient is 3/8-pi^2/2+6*zeta(3), and its complete first
moment is

    g1 = -13/4 + pi^2/2 - 2*zeta(3) = -0.719311605774...

This diagram class is not independently a conserving scalar probability
kernel. The source/pair-interference cancellation is discussed in
[Arbuzov, section 5.2](https://arxiv.org/html/hep-ph/0206036v2).
The code leaves the source intact and explicitly conditions the decay shape:

    K = G - g1*delta(1-z).

This follows by expanding E_raw(N)/E_raw(1) to the retained accuracy. The raw
source factor and removed first moment are logged. No weight is clipped or
fitted to PHOTOS/KKMC.

The formula, number normalization and numerical implementation are verified.
The outstanding physical question is whether this conditioned source sector,
together with our hard match, supplies exactly the photon-only logarithmic
coefficient relevant to the external generator observable. Normalization alone
cannot answer that question. Neither substituting a spacelike kernel nor
adding a pair-interference term merely to improve agreement would be justified.

### Explicit complete-radiator formulas

Set L=log(Q^2/m_mu^2), Lf=log(mu_frag^2/m_mu^2), ell=L-Lf, and let P_N, C_N,
K_N denote the plus moments of the distributions above. Then

    U_N = exp(2*a*Lf*P_N)
    R_LL  = U_N*(1+2*a*(C_N+ell*P_N))
    R_NLL = R_LL+2*a^2*Lf*U_N*K_N.

The complete first and second coefficients are

    r1 = 2*(L*P_N+C_N)
    r2_LL  = 2*Lf^2*P_N^2+4*Lf*P_N*(C_N+ell*P_N)
    r2_NLL = r2_LL+2*Lf*K_N.

At mu_frag=Q the retained coefficients of L^2 and L are 2*P_N^2 and
4*P_N*C_N+2*K_N. For N=2, P_N=-4/3, C_N=22/9 and K_N=3.51372243059795.
In particular r1_2=44/9-8*L/3 is an independently integrable first-order
identity. These equations specify what the current label means; they do not
include a complete constant r2.

The optional `quadratic_finite` mass-only prescription adds
2*a^2*U_N*(C_N+ell*P_N)^2. This is one finite O(alpha^2) completion, beyond
NLL. It preserves first order and the retained mass-log coefficients. It is
neither full NNLL nor a rigorous uncertainty bound.

## External radiation ingredients and their limits

### ReneSANCe/SANC: the sharp first-order check

The optional driver uses the pinned public ReneSANCe 1.3.11 SANC amplitudes:
`iqed=4`, `iew=0`, `iqcd=0`, with higher-order EW and QED evolution off and
LoopTools initialized. For the mass-only Z study, a massive angle-integrated
real-emission expression acts directly on F(real)-F(Born). The endpoint
subtraction supplies the unitary virtual term. Parent luminosities and mass
convolution use deterministic quadrature.

The upstream integrated NC routine does not honor `iqed` in the same way as
its differential counterpart. Our pinned source-preparation script selects
the eight FSR terms multiplying the current at s, preserving their coefficients;
it excludes the s-prime current and interference. The independent differential
MC check tests that selection. This is external matrix-element information,
not another use of our splitting functions. The convention follows the
[SANC–PHOTOS comparison](https://doi.org/10.1140/epjc/s10052-013-2625-1).

W uses massive real emission and virtual-plus-soft terms, with an independently
integrated inclusive normalization counterterm. The W+ neutrino/charged-lepton
ordering differs from W- and is handled explicitly. Soft slicing is varied.

The refined Z O(alpha) spectra differ by at most 1.16e-6 of the Born-bin rate
in the documented initial comparison. Its larger coarse-to-fine shifts were
1.41e-5 for fixed order and 1.15e-5 for NLL; that number is not a claim of
one-ppm final numerical certification. Direct differential versus integrated
Z MC has bin-pull RMS 1.013. W soft-cutoff pull RMS values are 1.101 and 0.987.
These are diagnostics, not global chi-square probabilities.

The separately labelled additive prediction

    SANC_NLO + (leptonff_resummed - leptonff_order_alpha)

uses the independent massive NLO anchor. It includes its finite-mass first-order
terms and our higher-order leading-power corrections. Agreement at NLO in that
prediction is imposed by matching; it cannot independently validate our
original first-order calculation. Those two comparisons remain separate.

### PHOTOS

The local reference is PHOTOS 3.56, exponentiation and interference on, pairs
off, fixed alpha, exponentiation tolerance 1e-7 and relative infrared cutoff
1e-5. Full production-dependent W/Z matrix-element weights are off because the
adapter supplies a standalone decay record. The legacy W correction is on.
This is the exact recorded configuration, rather than a generic statement
that “all MECs are off” or “PHOTOS is exact NLO.”

Finite hard radiation is also present in the analytic matched decay, so a
decay-level correction in a comparator is not automatically incompatible with
our formalism. The settings still differ in their detailed finite terms and
higher-order completion. A precision comparison requires documenting and
testing those differences, especially for W.

Four-million-decay cutoff-control samples per process reduce the cutoff by a
factor ten; the documented pull RMS range is 0.936–0.987. This finds no bias
resolved by those samples, rather than a zero-cutoff proof.

### PYTHIA

The reference is PYTHIA 8.318, a forced timelike decay shower with fixed alpha,
photon conversions and all other showers off. Decay matrix-element corrections
are on in the central comparison; MEC-off is a separate control. The lepton
shower cutoff is 1e-5 GeV, with dimensionful units unlike PHOTOS' relative
cutoff. Individual daughter shower scales are set explicitly. The adapter
follows the final muon through emissions with `iBotCopyId()` and checks that it
is final; a simple copy-chain lookup would return the wrong momentum.

The recorded cutoff pull RMS range is 1.022–1.034. MEC changes are measurable:
the Z fiducial correction in the earlier control is -1.154+/-0.017 percent
without MEC and -1.067+/-0.013 percent with MEC. These MC errors are not a
theory-uncertainty assignment. A shower is an independent algorithmic reference,
not automatically a fixed-order/NLL benchmark of identical content.

### KKMC

The fixed-Q adapter uses KKMCee v4.32.01, commit
`527ed9d6f8c0080a256e77ae0d7adbf4c5207eea`, without changing upstream
amplitudes. ISR, initial-final interference, weak corrections, real pairs,
the optional virtual-pair contribution and running alpha are off. Bare muon
momenta define the measurement. Each order's own inclusive weight normalizes
its moments and cumulatives.

EEX weights 72/73 and CEEX weights 251/252/253 have different residual-order
definitions. They are exponentiated predictions, not strict O(alpha) or
literal substitutes for our LL/NLL choices. The source's weight 74 does not
supply complete third-order FSR: its additional terms enhance ISR and do not
provide a full pure-three-photon FSR residual. Its apparent order label must
not be used to certify an alpha^3 FSR calculation. The formalism is described
in [the CEEX paper](https://arxiv.org/abs/hep-ph/0006359).

The defaults impose `MinMassCEEX=1 GeV`. The CEEX denominator excludes events
below that pair mass. Our leading-power radiator has wider support. The
first-order estimate of the omitted fraction is 0.70 ppm at the physical muon
mass, and 3.27 ppm at m_mu/100. The effect raises N>1 normalized moments,
opposite to the observed negative residual. This is not an exact all-order
threshold bound; identical support remains necessary for a definitive test.
The driver now exports the floor and pooling rejects different floors.

Photons below a CEEX technical threshold can be absent from its visible record
while their recoil remains in the muons. The adapter therefore bounds the
missing soft momentum rather than demanding exact closure of an incomplete
visible record. Radiation-off control uses EEX because this pinned CEEX release
does not support the all-radiation-off partition; CEEX has a near-zero-alpha
control instead.

### McMule: strict second-order shape

The new optional adapter uses McMule v0.8.0, commit
`de32e762ca44ab82345d4a511f7b64ab61ba9bd4`, with its pinned OpenLoops,
COLLIER and other submodules. It selects pure photonic muon-line real,
real–virtual and double-real pieces of a vector-current process. Vacuum
polarization, real pairs, ISR, IFI and weak corrections are excluded.
Finite lepton masses regulate collinear limits, and upstream FKS-l subtraction
handles soft cancellation; see the
[McMule formulation](https://arxiv.org/abs/2007.01654).

Writing rates as sigma=sigma0+alpha*sigma1+alpha^2*sigma2 and
M=sigma0+alpha*M1+alpha^2*M2, define I1=(M1-sigma1)/sigma0 and
I2=(M2-sigma2)/sigma0. The normalized result is

    R = 1+alpha*I1+alpha^2*(I2-(sigma1/sigma0)*I1).

For F(Born)=1, the no-real-photon double-virtual contribution cancels in I2.
The measurement F-1 must still act separately on every real and subtraction
counterevent. Integrating real–virtual and double-real independently without
their upstream subtraction, or forgetting the inclusive NLO counterterm,
would give the wrong answer. McMule sets internal alpha=1; the adapter restores
the relative powers of physical alpha after dividing by the Born rate.

The benchmark is a pure vector-current coefficient test. It is not a new full
Z NNLO electroweak calculation, nor an NNLO QCD run. It can test universal
mass-logarithm coefficients; finite-current and power effects must be kept
explicit when connecting it to Z phenomenology.

## Numerical implementation and error interpretation

The photonic backend uses log-loss cells, an endpoint atom and FFT transport.
The NLL insertion is signed. Probability clipping is not used. Number,
imaginary-part, omitted-tail and weighted-wrap diagnostics are distinct from
grid errors. The positive lower support of actual measurements matters when
undoing an energy tilt. Sharp bins require refinement even if sum rules and
wrap bounds pass.

The full-radiator audit uses independent continuum rational moments, an
independent dilogarithm/quadrature implementation, and direct double-plus
cumulative integrals. It checks 65,536 and 131,072 cells, several parent masses,
muon masses, couplings and fragmentation scales. Small-alpha extrapolation
extracts the **complete** second-order coefficient, with compensated summation
and a cubic-remainder check. A mass-log fit separately checks the implemented
L^2, L and constant pieces. There are 4,623 passing checks. This verifies the
displayed prescription, including its hard coefficient; it cannot by itself
prove an external-sector identification.

For common production or radiation streams, raw within-bin model covariance
is retained before forming ratios. Independent seeds are pooled before
normalization. Signed real/dipole production groups remain intact. Denominators
consistent with zero are withheld; unstable ratios are not converted into
apparently precise corrections. Deterministic quadrature nodes are not MC
events and do not produce statistical error bars.

The original comparison plot bands show **one-standard-deviation MC errors**.
The later leptonff accuracy plots separately display the three-point
fragmentation-scale envelope, with mu_frag/Q=1/2,1,2, and explicit numerical
refinement estimates. A scale envelope is not a statistical confidence interval
and does not cover all unknown terms. Grid/quadrature shifts are estimates,
not rigorous total numerical bounds. The finite-completion shift and generator
spread are separate diagnostics.

For KKMC, the full cross-order covariance is available within each observable.
The covariance between different moments or cumulative cuts is not exported;
their pulls cannot be summed into a global chi-square. McMule's native
histogram error is the scatter of iteration integrals. Independent seeds and
phase-space/subtraction controls therefore matter, particularly for heavy
tails. Its pieces use disjoint streams, and the normalization-product variance
is propagated. There is no invented cross-observable covariance.

## What the results establish

| Ingredient or claim | Evidence | Assessment and boundary |
|---|---|---|
| Whole-plus distribution conventions, endpoint actions and analytic first moments | Canonical analytic benchmark, C++ tests and independent Python integrals | Verified mathematical/implementation statements |
| Timelike kernel transcription and index convention | Source equation checks; independent moments; flavour and momentum identities | Strong implementation evidence; source class and physical-sector interpretation remain distinct |
| Leading-log evolution, one-insertion matching and complete expansion of our formula | Composition/expansion tests; independent continuum moments/cumulatives; 4,623 checks | Verified implementation of the stated prescription |
| Z normalized O(alpha) mass correction | Independent massive SANC real expression, differential MC closure, quadrature refinements | Strong external validation over the declared inclusive-angle Born-QCD setup |
| W real current and spin response | Ward/spin tests, massive anchors, 486 polarized RECOLA tree points (maximum reported difference 3.2e-14), SANC comparisons | Strong tree/current evidence; not precision certification of every fiducial W spectrum or full loop EW correction |
| Host grouping, tensor conventions, zero-alpha/native limits, archive replay | MCFM/DYTurbo native controls and covariance/replay tests | Working and tested interfaces; not universal proof for arbitrary host cuts/scales |
| Coupled NLL flavour evolution and thresholds | Conserving full-system moments, independent flavour references, matching/transport tests | Optional implementation is tested; it was not enabled in the photon-only comparison |
| Photonic NLL mass logs in the inclusive normalized vector moments | Independent McMule O(alpha^2) scan, Q=10–1000 GeV; constant-remainder fits and FKS/technical controls | Supported at the measured precision; not a certification of all fiducial observables or all-order completions |
| Total phenomenological uncertainty and tenfold smaller errors in every fine bin | Partial scale, seed, grid, domain and cutoff studies | Not established for all requested observables, especially W |

The installed library has C++, C, Fortran and Python interfaces, grouped
production archives, replay, configurable measurement bins, cached responses
and parallel evaluation. These improve reuse but do not create an accurate
production sample or recover FSR migration from a histogram that has already
lost the needed kinematics.

## Resolved effects versus unresolved interpretation

In the leptonff-only Born-QCD Z sample, the integrated 90–92 GeV bin has:

| Prediction | Rate [pb] | Correction to common Born |
|---|---:|---:|
| Born | 677.934753 | — |
| O(alpha) | 583.174642 | -13.9778% |
| Matched photonic LL | 588.209864 | -13.2350% |
| Matched photonic NLL | 588.291115 | -13.2230% |

LL changes the first-order result by +0.8634% relative; NLL changes LL by
+0.013813% relative. Small visual separation of LL and NLL on a broad plot is
therefore expected. It is not evidence that the NLL code is inactive. The
refined NLL-minus-LL bin differences are resolved numerically in all 264 mass
bins; the largest recorded numerical estimate for that ratio is 3.34e-9.
That particularly small paired numerical effect is not the precision of the
absolute physical cross section. At first order, scale cancellation is checked
to about 3.3e-15 in the recorded ratios.

At fixed Q=M_Z, KKMC CEEX2 minus leptonff NLL gives for N=2
-86.06+/-7.13 ppm. At half alpha it is -20.30+/-4.90 ppm; at twice alpha,
-441.11+/-10.71 ppm. Reducing the muon mass by a factor 100 gives
-326.38+/-9.98 ppm at physical alpha. The additional targeted campaign contains
850 million trials, beyond the existing 200-million-trial physical points and
100-million-trial soft-cutoff control. The tenfold cutoff change has maximum
pull about 1.02 in the recorded observables/orders.

The mass dependence rules out the tested single finite-completion ansatz as
the sole explanation of the full resummed difference. It does not uniquely
prove a missing alpha^2*L coefficient: coupling fits admitting higher powers
become poorly constrained, and different higher-order completions can acquire
mass dependence at finite alpha. A successful internal coefficient check is
also not a substitute for this physical test. The
[residual investigation](z_radiator_residuals.md) records the actual fitted
coefficients and their truncation dependence.

The new strict fixed-order check changes that assessment in a useful way.
At M_Z the N=2 alpha^2 contribution from McMule minus our NLL expansion is
-28.073+/-0.719 ppm. Across six Q values, a constant fit gives
-28.660+/-0.317 ppm with chi-square 3.17 for five degrees of freedom. The
corresponding constants for N=3,5,10 are -17.526+/-0.432,
+5.558+/-0.604 and +47.551+/-1.020 ppm. Additional L/L^2 terms are not resolved;
the analysis retains their correlated fit errors and lower-Q power sensitivity.
Subtraction and tighter technical controls have a largest pull of 1.56.

This supplies independent support for the **mass-logarithm coefficients** in
the tested observable and identifies a finite second-order shape remainder.
It does not identify every term in the resummed KKMC residual. The largest-Q
points are less precise, and near-endpoint cumulative cuts have much larger
errors than these moments. No fitted coefficient has been inserted into the
production prediction. Full finite NNLO matching remains a separate task.

## Uncertainties and ingredients still missing

The [research extension notes](manual/extensions.md) develop possible NNLO
matching, endpoint and soft/collinear improvements, including joint cone/veto
resummation for dressed electrons. They distinguish universal fragmentation
from process-dependent Z/W hard terms and retain the current accuracy claims.

1. **Perturbative/sector closure.** The strict fixed-order moments now support
   our NLL mass-log content within their fit uncertainties. Extend that evidence
   to the full desired observables, align CEEX support exactly, and interpret
   the remaining resummed differences and finite alpha^2 function. The present
   fixed-order check is not a full all-order or fiducial certification. Do not
   tune the published kernel to a shower.
2. **Higher QED orders.** Complete two-loop finite hard/boundary terms and full
   NNLL resummation are absent. Their uncertainty is not determined by one
   scale envelope or one quadratic finite completion. Soft/collinear endpoint
   and very small-z regions need observable-specific scrutiny.
3. **Pairs and running for physical precision.** The photon-only comparison
   deliberately excludes real pairs and fermion-loop running. The library
   has separate coupled lepton/flavour/threshold prescriptions. A physical
   prediction needs a specified extra-pair selection and electromagnetic
   input scheme; a photon-only generator spread does not bound those effects.
4. **Hadronic content.** HVP inputs and a pointlike quark prescription exist as
   optional extensions. HVP does not determine real hadronic-pair
   fragmentation. A fitted hadronic boundary and complete mixed QCD/QED
   fragmentation evolution are not supplied. They are not prerequisites for
   the restricted photon-only benchmark.
5. **Production and electroweak content.** The current precision Z comparison
   has Born QCD production. PDF/alpha_s variations, QCD scale variations,
   recoil/resummation choices and matching to the chosen experimental setup
   are not included in its QED error bands. Full weak corrections, QED ISR,
   IFI, photon-initiated production and complete mixed QCD–EW corrections are
   outside this decay-only calculation. A factorized QCD production times QED
   decay prediction contains selected mixed products, not the complete mixed
   perturbative correction.
6. **Finite mass and recoil.** The analytic decay is leading power in m_mu/Q,
   with explicit massive external references. Its collinear momentum map is
   not a fully exclusive multigamma event generator. Generic fiducial/recoil
   observables, pair threshold, asymmetric/angular/parity-odd measurements
   and photon observables need separate contracts and validation. The current
   Z hard response has a parity-even charge-symmetric scope.
7. **Precision W phenomenology.** NLO+NLL QCD W coupling works in the DYTurbo
   interface, but a full high-statistics fine-bin W comparison, hard-remainder
   convergence and a controlled uncertainty budget are still needed. Born
   baselines are useful validation; they do not describe QCD recoil near a
   Jacobian edge at precision level.
8. **Observable/data definition.** Bare and post-FSR often coincide only after
   specifying promptness, extra pairs, stable-particle muons, photon treatment,
   parent domain and fiducial cuts. No selected measurement's unfolding
   definition has been certified here. No detector resolution is modeled.

Numerical errors must also be reduced where needed: independent production
statistics, hard/angular sampling, integration iterations, FF lattice,
quadrature, omitted domain and technical cutoffs. More statistics cannot
repair a wrong normalization, a mismatched sector or a missing perturbative
coefficient. Conversely, unresolved MC differences should not be interpreted
as theory discrepancies.

## Useful additional validation observables

Normalized fixed-Q moments and cumulative cuts remain especially useful: they
remove PDFs and rate conventions, admit analytic expansions and probe
different radiation regions. Inclusive no-cut rate conservation is necessary
but weak, because our convention enforces it. A massive first-order shape and
an independent fixed-order second coefficient are stronger tests.

Dressed-muon cone scans can test the approach toward collinear safety, but
they need photon-angle/recombination information beyond a purely collinear
fragmentation action; they are not automatically predictions of the present
mass-only reduction. Bare-to-Born migration, cumulative energy loss and
matched fixed-parent samples are available without pretending to generate
exclusive photons. Z pT/eta cuts, W pT/eta/mT and different parent/recoil bins
are subsequent tests once their numerical and angular contracts are met.

## Reproducibility and current evidence files

The main evidence is in [the first-order benchmark](born_fsr_reference.md),
[the accuracy-scale snapshot](../benchmarks/z_accuracy_scales.json),
[the full-radiator audit](z_radiator_validation.md),
[the residual snapshot](../benchmarks/z_radiator_residuals.json),
and [the new fixed-order guide](mcmule_fixed_order.md).
The source and exact commands for each external reference live in
`adapters/renesance`, `adapters/kkmcee` and `adapters/mcmule`. The production
interfaces are documented in [the DYTurbo guide](dyturbo_integration.md),
[the W coupling guide](dyturbo_w_decay.md) and [the manual](manual/index.md).

External sources, PDF grids and linked executables retain their own licence
and citation requirements. They are not redistributed under the library's
MIT licence. The McMule client retains a documented explicit two-body dispatch
workaround; no amplitude or FKS formula is edited. Its importance maps are
independently checked against native phase-space volumes and kinematics.

The obsolete cached Python 3.14 path that prevented a previous CTest invocation
has been replaced by a fresh validation environment and a reconfigured build.
Fresh test outcomes and the final fixed-order assessment are recorded in the
new benchmark guide, rather than retroactively changing historical logs.
