# Photonic hard dressing on the W/Z NLO QCD host

Version 0.24 also supplies a separate fixed-alpha, unitary-conditioned
[photonic NLL mode](photonic_nll.md), with a single time-like NLO-kernel
insertion. The LL branch described below and its old results are preserved.

Version 0.17 dresses the spin-correlated W and charge-symmetric Z finite decay
term with pure photonic LL evolution. It also evolves the MSbar boundary,
retaining **one finite insertion globally**. The leading-power unitary decay,
spin tomography and common real/dipole random streams retain the
[v0.16 contract](spin_decay_matching.md).

This completes a missing ingredient for NLL, not a complete physical NLL
prediction: this branch still has no P1, running, pairs or physical flavour
thresholds. The coupled one-family NLL kernel reference is a separate branch.

## Operator and overlap

Write a=alpha/(2*pi), L=log(mu_frag^2/mu0^2), U=exp(a L P0), and B=U*d1.
For n charged legs the Born-parent action is

```
FF = U^n + a sum_r B_r product_(s!=r) U_s
     + a*(D1_continuum - C1_FF_lattice).
```

C1_FF_lattice includes both the actual LL generator and the original lattice
d1. For Z, BxB is absent. The explicit C2 is also supplied: its single-leg
coefficient is L^2 P0*P0/2 + L P0*d1, and its two-leg mixed term is
C1xC1-d1xd1. No constant d1xd1 term is silently promoted to this accuracy.

The dressed hard coefficient is constructed as

```
H_dressed = H_bare_continuum
          + H_projected_lattice[U^n] - H_projected_lattice[delta^n]
          + G[U^n] - G[delta^n]
R_dressed = FF + a*H_dressed.
```

The bare coefficient is the v0.16 continuum projected term plus its finite
geometric integral at the declared numerical cutoffs. The subtraction uses
exactly the projected lattice that was evolved. It does not substitute a
continuum approximation for a different numerical overlap. Hard dressing
and boundary dressing therefore start beyond O(alpha), even on a coarse
grid. With identity evolution, both hard dressing shifts vanish exactly and
the result reduces to the previous NLO operator. With zero coupling, all
radiative columns vanish and both W and Z reduce to the Born measurement.

Only U acts on hard-event charged leptons inside H1; using U+a*B there would
insert H1*d1. Photons are unobserved and do not convert to muons in this
photonic branch. For W the neutrino is not fragmented and visible collinear
remnants are not assigned to missing momentum.

## Projection and real kinematics

`PhotonicLabDressing` projects the per-leg vector hard kernel at M=mu_frag,
the additional charged-W term, and P0 onto the common log-loss lattice.
Each whole-plus endpoint is retained, with no clipping or renormalisation.
The hard source is transported by zero-padded causal convolution: a high-loss
cell cannot wrap into the physical endpoint. Existing LL/boundary periodic
wrap and tail controls still apply to their respective evolutions.

For off-shell Z/gamma* kinematics the hard kernel is reconstructed as
`h(M=mu_frag)+log(M^2/mu_frag^2)*P0`. M is the original hard parent mass,
not the mass after fragmenting a daughter. Both Z legs are evolved in the
hard projection. All laboratory pT, eta and invariant-mass bounds are applied
to the resulting signed tensor measure before multiplication by the host
weight.

The geometric operator evaluates its bare and dressing-shift tensors from
one current, rotation and local collinear projection. Its evolved measurement
uses exact signed LL prefix sums for W and the two-leg prefix contraction for
Z. The same evolved measurement is used for the real term and every local
subtraction. This avoids independent fragmentation draws in cancelling
terms. At identity evolution the geometric difference is zero point by point.

## Host columns and statistics

`--dress-hard` requires `--photonic-ll --hard-samples N`. Earlier columns
retain their definitions and reproduce the v0.16 nominal results bitwise.
The new rates are:

- `boundary_dressing`: boundary-evolved FF minus the earlier FF baseline.
- `projected_dressing`, `geometry_dressing`: a times the hard dressing shifts.
- `dressing_shift`: the sum of those three contributions.
- `dressed_matched`: previous additive result plus dressing_shift.
- `photos_minus_dressed`: a directly accumulated paired comparison.
- `dressed_minus_order1`: the directly accumulated higher-order remainder.

The optional raw iteration audit includes all 26 rate columns. Common
iteration means retain algebraic closure and a rate covariance estimate;
independently averaged native component histograms must not be subtracted to
construct an error. Ten iteration chunks per QCD contribution remain a
limited covariance diagnostic, not a heavy-tail or cross-bin precision bound.

See the [benchmark](../benchmarks/photonic_dressing_reference.md) for runs,
controls, actual shifts and provenance. Physical multi-flavour/coupling/pair
content, all hard-parent channels in the coupled NLL branch, P1/running in a
consistently defined photonic comparison, differential hard histograms and a
controlled phenomenology error budget remain separate work.
