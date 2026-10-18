# Spin-correlated finite decay matching

Version 0.16 supplies leading-power, unitarily normalised O(alpha) W and
neutral-vector decay operators, and an optional additive MCFM rate pilot.
The W current includes photon emission from the charged W as well as the
charged lepton. The neutral current includes both charged legs. This is a
factorised decay correction on NLO QCD production, not complete NLO EW or
complete mixed QCD–EW corrections.

## Normalisation and finite subtraction

Write `a=alpha/(2*pi)` and `P(z)=(1+z*z)/(1-z)`. Here `z` is the invariant
mass squared of the two final fermions divided by the decaying boson mass
squared. In the collinear limit it is also the radiating lepton fraction.
All distributions below use a whole-plus prescription.

The independently integrated real currents have the leading-power spectra

```
r_W(z) = P(z)*(log(M*M/m_mu/m_mu)+log(z)-1)
         - 2*(1+z+z*z)/(3*(1-z))
r_V(z) = 2*P(z)*(log(M*M/m_mu/m_mu)+log(z)-1).
```

The W spectrum is also checked against the analytic public SANC decay
spectrum recorded in the neighbouring `w-weak-reweight/TECHNICAL_NOTE.md`.
The real W graph follows the ordinary WWgamma vertex and unitary-gauge
propagator; the photon Ward cancellation is checked before squaring. This is
not a deletion of photon diagrams from a full off-shell production amplitude.
For context on gauge-consistent pole corrections see
[Dittmaier–Krämer](https://arxiv.org/abs/hep-ph/0109062), and on real/virtual
photon corrections to decays see the
[Sherpa decay calculation](https://link.springer.com/article/10.1140/epjc/s10052-019-6614-x).

The inclusive decay correction is defined to vanish: the finite inclusive
normalisation belongs to the chosen branching ratio. In particular this
removes the `3/2` vector-current inclusive coefficient used by the older
unnormalised vector reference. For W it corresponds to the FSR-normalised
convention, rather than implicitly using the public SANC muPW-dependent
inclusive remainder. This convention must also be used by any later weak or
branching-ratio correction. PHOTOS is used as a unitary decay generator.

Subtracting the implemented MSbar D1 gives the per-leg projected hard kernels

```
h_V(z) = P(z)*(log(M*M/mu_frag/mu_frag)+log(z)+2*log(1-z))
h_W(z) = h_V(z) - 2*(1+z+z*z)/(3*(1-z)).
```

The whole-plus endpoint, including the finite W term, is retained. Rate
windows are integrated before the host applies its event weight; no fitted
PHOTOS correction enters H1. The W second invariant-mass moment is
`11/3-4*log(M*M/m_mu/m_mu)/3`; the neutral one is
`44/9-8*log(M*M/m_mu/m_mu)/3`. Hard+FF scale cancellation against these
independent moments is tested for mu_frag=M/2,M,2M.

## Finite angular remainder

`spin_radiative_point` returns a Cartesian 3x3 Hermitian tensor R in the boson
rest frame, with measure `dz dc dR`, where `dR` is normalised SO(3) Haar
measure. `spin_born_tensor` returns B(n), normalised so `<B>_sphere=I` and
`Tr B=3`. A production tensor rho contracts as `Tr(rho B)` or `Tr(rho R)`.
There is no boson spin average hidden in a polarised call.

Set `r=Tr(R)/3`. For W project the charged-lepton+photon cluster onto a
massless two-body Born direction opposite the neutrino. For the neutral
current partition the two collinear regions with `(1+c)/2` and `(1-c)/2`,
and use the corresponding emitting cluster. Each neutral sector integrates
to half the real energy spectrum. The geometric operator is

```
G[F] = R*F(real) - r*sum_s w_s B(n_s)*F(collinear_s(z)).
```

Its soft and collinear leading terms cancel locally. Integration over global
orientation removes G for any measurement depending only on the two-fermion
mass. Thus `H1[F] = projected_hard[F] + integral G[F]` has the required
spin information and finite noncollinear migrations. Numerical soft,
collinear and very-small-z exclusions are refinement parameters on G, not
factorisation counterterms. Their residual shifts must be measured.

An optional orientation-independent reference subtracts
`(R-r*sum_s w_s B(n_s))*F_ref`, whose orientation integral is exactly zero.
The host uses 1 for W and the invariant pair-mass selection for Z. This makes
fully inclusive angular measurements cancel point by point and reduces
variance. The reference must be independent of the sampled orientation;
using a direction-dependent Born cut here would change the answer.

## Extracting spin from a scalar NLO host

The optional host correction uses angular tomography rather than dividing by
a Born weight. For W, with its charge-dependent sign s,

```
B(n) = 3/2*(I-n*n^T+s*i*epsilon(n))
K(n) = 2*I-5*n*n^T+s*i*epsilon(n)
<Tr(rho B(n))*K(n)>_sphere = rho.
```

The signed K has trace 1; it need not be positive. The host samples an
independent radiative orientation, then contracts G with K of the original
charged-lepton direction. After integration over the original Born decay
this recovers the full production spin tensor, including QCD recoil and
negative weights. Cubature checks reconstruct all nine Hermitian components
for both W charges. This identity requires full Born angular coverage at
fixed production kinematics: host lepton cuts must be disabled, and
production selections/scales must not select the original decay direction.
The provided runner uses those conditions. It is a rate estimator; one cannot
label the hard correction with the original Born eta and obtain a physical
radiative eta spectrum. New hard columns therefore have rate histograms only.

For gamma*/Z the present adapter reconstructs only the real symmetric tensor,
using `K_even=2I-5*n*n^T`. Massless vector and axial currents have the same
parity-even decay tensor. The parity-odd part vanishes for the charge-symmetric
pair cuts in this pilot, including off-shell photon/Z interference. This does
not implement a forward-backward asymmetry or a charge-asymmetric Z spectrum.

Real and mapped QCD subtraction terms use the same independent photon random
stream, seeded from the original real lepton momenta. Counterevents must never
receive independent radiation samples. Hermitian pairs are contracted together
to avoid spurious imaginary residuals from soft cancellation in double
precision. Host-level subtraction/refinement checks remain necessary in
addition to the algebraic identities.

## Fixed-alpha photonic comparison

`LEPTONFF_PHOTONIC_LL=1` selects fixed-alpha non-singlet LL evolution from a
unit endpoint, with no photon-to-pair feed-back, wrong-sign pairs, running or
P1 kernel. Its actual LL lattice coefficient is subtracted and the same
continuum MSbar D1 is restored. The additive prediction is

```
FF_continuum_matched = U_LL + a*(D1_continuum - U1_lattice)
decay_order1         = Born + a*(D1_continuum + H1)
decay_matched        = FF_continuum_matched + a*H1.
```

The numerical alpha is the same alpha(0) in this fixed photonic branch and
PHOTOS. Calling this a finite conversion of the pair-inclusive MSbar NLL
calculation would be incorrect; the distinct photonic branch avoids that
comparison. The default coupled one-family FF mode is retained.

PHOTOS first-order mode disables exponentiation and multiple emission. At
physical alpha and cutoff 1e-5 its *crude proposal* exceeds the library's bound
in part of the high-mass Z generation tail; even 1e-3 is insufficient across
the complete tail. No probability is clipped and those attempts are retained
as failed runs. The supported benchmark instead samples at alpha/s, s=4,
and returns `F_B+s*(F_after-F_B)`. With W ME on / Z ME off the single-photon
transition is linear in alpha, so this is an unbiased signed estimator of the
requested first-order operator. Events with more than one photon are rejected.
The physical alpha and 1e-5 cutoff are unchanged. A paired statistical closure
against direct physical-alpha sampling at a supported rest-frame point passes
for W+, W- and Z (50,000 events per branch and process). See
`adapters/mcfm/check_photos_linearity.py`.

The proposal scale is explicit in each manifest and must equal one for
exponentiation. Exponentiated runs use physical alpha directly and the same
1e-5 cutoff. The local setup has W ME correction on, Z ME correction off,
and interference on; agreement with its angular approximation is a measured
benchmark, not an assumption of identical amplitudes.

The additive result is NLO decay + LL, not full NLL. Evolving H1 on all hard
parents, boundary-generated subleading terms in this photonic branch,
physical multi-flavour/pair content, differential hard histograms, wider phase-space validation and a controlled
statistical/numerical error budget remain required for a production NLL phenomenology prediction.


## Independent polarised amplitude check

The optional RECOLA 1.5.1 driver in `adapters/recola/` evaluates all three
boson helicities at 486 rotated/boosted real-emission points. Both W charges
are compared directly; the physical Z amplitude is symmetrised in the two
lepton momenta to isolate the parity-even current used by this adapter.
External leptons are massless, widths are zero, and the reference normalises
to each process's unpolarised Born amplitude with the explicit phase-space
factor `M*M*(1-z)/(4*e*e)`. No fitted normalisation or diagram subtraction is
used. Maximum relative differences are 3.16e-14 (W-), 3.13e-14 (W+) and
5.59e-15 (Z). The complete grid, commands and library/executable hashes are
in `results/recola_spin_checked/`. This independently checks the real
polarisation structure; the unitary inclusive virtual normalisation remains
the explicitly chosen decay convention.


## Rate-estimator audit

MCFM normally combines iterations with a different inverse-variance weight
for each histogram. These weights are estimated from the same samples as the
means; small finite hard components and cut shoulders can therefore have
slow or biased finite-sample convergence. Component rows do not inherit the
paired errors of their sum. The W- projected hard component triggered a
3.16-sigma subtraction-parameter diagnostic in the initial pilot, although
the complete matched rate passed. An independent seed with four times the
calls reduced that component diagnostic to 1.90 sigma. This does not erase
the initial flag or certify its component precision.

`build_adapter.py --iteration-audit` optionally instruments an isolated copy
of MCFMStorage immediately before its raw iteration histograms are reset.
It logs stage, QCD part, phase-space channel, iteration, calls, rate and
within-iteration variance for all rate columns, retaining zeros and the
original correlated real/dipole sums. The original installation and its
integration/weighting are untouched. `analyse_iterations.py` reconstructs
the original weighted estimates, then forms common equal-iteration means
within fixed-call contributions and checks algebraic closure. A full rate
covariance is estimated from the iteration means; with only ten chunks per
contribution it is a diagnostic rather than a precision covariance matrix.
No cross-bin covariance or heavy-tail error bound is inferred.


The common-iteration W- subtraction checks give pulls -0.66 (projected),
+0.41 (geometry), +0.24 (hard sum), -0.18 (matched) and -0.66 (paired PHOTOS
difference). Component closure is within 1.2e-10 pb; the original MCFM
weighted results are reconstructed to printed precision. The common-iteration
errors are larger and retain the rare-event contributions. Passing these
checks establishes a useful pilot, not a certified precision error budget.
