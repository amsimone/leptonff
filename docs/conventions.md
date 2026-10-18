# Implementation conventions

All energies and masses are in GeV. Four-vectors use `(E, px, py, pz)` and metric
`(+---)`. The collinear FF fixture uses massless hard/observed four-vectors: its
muon mass enters the perturbative FF logarithms. Separate massive real/virtual
references retain mass in the currents and kinematics. Collinear rescaling is
`p = z*k`, `0 ≤ z ≤ 1`.
Pseudorapidity is undefined for a zero three-momentum.

The FF interface is `parent -> observed`. The paper's `D[H/a]` maps to
`D[a -> H]`. Hard phase-space weights multiply the measurement evaluated at
rescaled momenta; no extra `1/z` belongs in this event-functional convention.
The current basis is exactly `{mu-, mu+, gamma}` with unit lepton charge.

Write `D = I + a D1 + O(a²)`, `a = alpha/(2*pi)`. Coefficient functions do not
include `a`; `fragmentation_order1` does. `initial_condition_coefficient` is
evaluated at `mu_0`. Scale transport is explicitly truncated:

```text
D1(mu_frag) = d1(mu_0) + log(mu_frag²/mu_0²) P0
```

This scale-transport API is fixed-coupling O(α), distinct from the all-order
non-singlet LL API in [ll_evolution.md](ll_evolution.md). `mu_frag`
and `mu_0` are distinct from future QCD `mu_R` and proton `mu_pdf` scales.
MS̄ is the only implemented FF scheme. The fixture has no EW Born normalisation,
PDFs or flavour thresholds. This describes the original fixed-order fixture;
the separate [coupled solver](moment_evolution.md) supplies MSbar running, and
the [hard-matched moments](photon_pair_matching.md) use the same MSbar input for
FFs and hard alpha(M). Neither supplies physical multi-flavour or exclusive
extra-pair selection.

## Distribution action

The plus prescription acts on the entire kernel, including a z-dependent
numerator. For `0 ≤ c < 1`:

```text
∫[c,1] [g(z)]+ f(z) dz
 = ∫[c,1] g(z)(f(z)-f(1)) dz - f(1) ∫[0,c] g(z) dz.
```

Measurements must have the endpoint limit represented by `f(1)` and sufficient
regularity for this subtracted integral to exist. Endpoint-only or discontinuous
endpoint measurements need a separate definition. A zero-width `[1,1]` action
is rejected. At nonzero lower bounds the endpoint delta is included.

`apply_pair` nests the two distribution actions; this includes all cross terms
and the double-endpoint term. For cuts, supply outer breakpoints and inner
breakpoints as a function of the outer z. Merely using an adaptive rule cannot
guarantee detection of a narrow bin. At fixed O(α), use `apply_pair_order1`
instead of multiplying two fully truncated `I+aD1` distributions.

Real Mellin moments mean `∫ dz z^(N-1) D(z)`. No complex moments are implemented.
The lepton-to-photon channel requires N > 1 because of its small-z singularity;
its photon number moment is deliberately rejected. Other O(alpha) channels
accept N > 0. Full NLO singlet/photon channels can also have small-z poles;
the complete matrix solver requires N>1. The analytic valence combination
supports its conserving N=1 check. Sum the N=2 moment over all daughters for
a fixed parent.

## Numerical controls and reference geometry

The standalone integrator is adaptive embedded Gauss–Kronrod 7/15 with explicit
breakpoints. It refines the interval with the largest estimated error until the
sum of estimates meets `max(abs_tol, rel_tol*abs(value))`. Non-finite samples,
resolution exhaustion and interval-budget exhaustion fail explicitly. These
estimates are diagnostics, not rigorous bounds. Distribution actions can combine
multiple integrals; nested pair actions do not propagate a global error bound.
Refine tolerances and compare with independent references for reported results.

`Real` is `long double`, but some platforms, including Apple ARM, implement it
with 53 mantissa bits. The driver records the actual precision.

The parent-at-rest isotropic fixture is normalised to one inclusive decay.
With `r = 2*pT_min/(M*min(z1,z2))` and `m = M*sqrt(z1*z2)`, the fiducial
acceptance is zero outside the strict mass window or for `r ≥ 1`; otherwise it is
`min(tanh(eta_max), sqrt(1-r²))`. This angular formula provides an exact
measurement for the longitudinal FF convolution. It has no DY spin information.

The fixed-order driver reports `Born + a*(D1⊗I + I⊗D1)[acceptance]`. Missing hard decay
terms mean `mu_frag` dependence is not cancelled in this reference. Starting-scale
dependence does cancel through O(α) because the boundary logarithms are retained.
The refinement shift is a numerical check, not a theory uncertainty band.
