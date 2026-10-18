# Fiducial vector decay and additive LL matching

Version 0.5 extends the unpolarised vector current at rest to the symmetric
`FiducialCuts` pT, pseudorapidity and bare-pair mass cuts. It supplies the finite
hard coefficient in MSbar, a massive fixed-order reference, and additive
matching to the existing LL plus single-boundary result. Normalisation remains
massless total vector Born. This is a restricted photonic decay calculation;
it does not provide a DY spin density, boosted parent, axial mass terms, running
coupling, pair channels or an NLL evolution kernel.

## Angular measurement

For each three-body configuration, rotational invariance makes the overall
orientation uniform. A cut on either muon becomes `|cos(theta_i)| < b_i`, where

```text
b_i = min(tanh(eta_max), sqrt(1-(pt_min/|p_i|)^2)).
```

The acceptance is the area shared by two spherical bands, divided by 4*pi.
The implementation evaluates their complementary spherical-cap intersections
analytically. It uses the vector cross product to determine the opening angle,
which retains precision for a nearly back-to-back pair. Tests compare this
formula with a distinct integration over the allowed azimuth interval. A
vanishing muon momentum or failure of the pT/mass cuts gives zero acceptance.

Only symmetric muon cuts and the unpolarised parent-at-rest model are supported.
The analytic averaging is an observable-specific projection, rather than an
export of fully differential weighted hard events.

## Local finite subtraction

Let z=m_mumu^2/M^2, c be the mu-minus pair-rest-frame cosine to the photon,
A(z,c) the rotational acceptance of the real event, and B(z) the isotropic
two-body acceptance after rescaling one leg by z. In either collinear limit
A(z,+/-1)=B(z). Symmetric cuts make the two limits equivalent. The massless
real density, relative to a=alpha/(2*pi), is

```text
R0(z,c) = ((1+z)^2 + (1-z)^2*c^2) / ((1-z)*(1-c^2)).
G[cuts] = integral dz dc R0(z,c) * (A(z,c)-B(z)).
H1[cuts] = H1_mass_spectrum[B] + G[cuts].
C1[cuts] = H1[cuts] + 2*D1[B].
```

The hard mass-spectrum distribution is defined in
[the previous matching note](vector_mass_matching.md). Its plus prescription
includes the finite virtual term. The difference defining G cancels collinear
singularities locally; its soft limit is integrable. The FF term equals the
existing `decay_ff_coefficient_order1`. Thus the hard and FF fragmentation-scale
logarithms cancel through O(alpha) for the fiducial observable too.

The scalar massive density follows the energy-fraction trace in
[Bystritskiy et al., hep-ph/0505236v2](https://arxiv.org/html/hep-ph/0505236v2),
eq. (3), retaining its mass-dependent double pole. An independent implementation
with explicit Dirac matrices checks it pointwise, including the massless limit.
Dropping the double pole would lose a finite collinear contribution.

## Massive reference

Massive Born and the integrated virtual correction are rotationally invariant.
The massive reference therefore combines the exact inclusive coefficient from
the earlier module with a finite real difference. With Bm(z)=Am(z,1) and
ABm the massive Born angular acceptance,

```text
Cm[cuts] = ABm * Cm_inclusive
           + integral dz Sm(z) * (Bm(z)-ABm)
           + integral dz dc Rm(z,c) * (Am(z,c)-Bm(z)).
```

Sm is the exact angle-integrated real spectrum. The rate is the massive
fiducial Born rate plus a*Cm. The angular integration uses a hyperbolic map to
resolve both mass-regulated collinear peaks. This reference makes a genuine
mass scan possible with the same cuts. Removing angular cuts reproduces the
mass-only module for inclusive and real-tail windows.

## Numerical regulators and matching

G omits z>1-epsilon and |c|>1-delta. These are numerical slivers of a finite
integral, not factorisation regulators. The massive calculation omits the
corresponding finite soft sliver but needs no angular cutoff. The C++ API
exposes these parameters and the inner/outer quadrature controls through
`HardDecayOptions`. The driver reports separate soft, collinear and quadrature
refinement shifts. Its defaults are tested fixtures, not error guarantees
for arbitrary narrow cuts. Tied eta-only thresholds give a collinear-sliver
error proportional to sqrt(delta) and can have linear mass corrections.
Their separate C++ convergence fixture uses delta=1e-10 and epsilon=1e-7;
the pT fixture's faster convergence must not be assumed. The strict upper
mass cut at M is rejected because
it removes Born but accepts arbitrarily soft radiation.

The combined prescription is

```text
matched = R_FF_lattice
          + a*(C1_FF_continuum - C1_FF_lattice)
          + a*H1_fiducial.
```

The existing FF matcher already implements the first two terms, including the
complete lattice overlap and at most one boundary insertion across both legs.
The new matcher adds H1 exactly once and rejects a cached hard coefficient whose
mass, cuts or QED parameters differ from the evolution configuration. Coupling
scans verify that matched minus full fixed order starts at alpha squared. With
no evolution it equals full fixed order; zero coupling gives Born.

This is additive NLO+LL for the stated decay model, with selected boundary
terms beyond LL. H1 is **not dressed with LL evolution**: hard-generated
alpha^n L^(n-1) terms are absent. It is not NLL and should not be presented as
the plan's final H1-convolved resummation. Version 0.6 supplies the separate
[LL hard-dressing extension](hard_dressing.md). P1/running/mixing remain
separate work. The angular and massive checks advance G2 within the
reference vector model; the physical spin-correlated Z gate remains open.
