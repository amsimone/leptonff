# Coupled spin-correlated W/Z matching

Version 0.18 joins the physical W± and parity-even, charge-symmetric Z/gamma*
finite decay operators to the existing **one-family MSbar NLL** evolution.
The kernels P0/P1, two-loop running, d1 boundaries, charge-flipped primaries
and photon-parent conversion all participate. This is a factorised decay on
an NLO QCD host. Physical multi-flavour thresholds, a finite conversion to an
on-shell input scheme, complete EW/mixed corrections and a phenomenology
error budget remain separate work.

## Measurement contract

W counts inclusive muons of the primary W decay charge. The neutrino is
unchanged by fragmentation; visible remnants do not enter missing momentum.
Z counts inclusive opposite-charge pairs from distinct hard directions with
positive observed pair mass, symmetric pT/eta cuts and no extra-pair veto.
These are identified inclusive densities once conversion is enabled, not
necessarily probabilities or an exclusive exactly-two-muon event selection.

Write S=D_mu->same, O=D_mu->opposite, G=D_gamma->mu_minus=D_gamma->mu_plus.
On the three-body hard event the measurements are

```
W: F(mu,nu; S) + F(gamma,nu; G)
Z: F(mu-,mu+; S,S) + F(mu-,mu+; O,O)
   + sum_(i=mu-,mu+) [F(gamma,i; G,S) + F(gamma,i; G,O)].
```

A pair from one unresolved parent direction has zero invariant mass and is
excluded by the positive pair-mass cut. No fragmentation fraction is drawn:
signed prefix contractions evaluate every term on both the real event and
its collinear projection. The energy tilt is undone only above positive
measurement support; no number normalisation or weight clipping is applied.
The neutral operator remains restricted to charge-symmetric measurements.

## Photon subtraction and projection

Retain the [spin operator](spin_decay_matching.md)'s z=m_fermion_pair^2/M^2
cluster projection and its two neutral sectors. In a collinear sector, the
emitting muon carries z and the photon x=1-z; the other Born daughter is
unchanged. For the photon channel the leading-power real spectrum is the
corresponding per-leg r(1-x). Subtract the same MSbar boundary/transport used
by the evolution,

```
D1_mu->gamma(x) = (1+(1-x)^2)/x * (log(mu_frag^2/m_mu^2)-2 log(x)-1).
h_gamma(x) = h_mu(1-x), REGULAR, with no plus prescription or delta.
```

Here h_mu is the full per-leg hard kernel, including the finite W term.
This photon kernel has a divergent number integral at x=0; its energy moment
is finite. Per charged direction, the unitary hard source obeys
`integral z*h_mu(z) + integral x*h_gamma(x) = 0`.

The geometric tensor is the same real R minus the sector-weighted Born
projection, now with all three momenta passed to the measurement. It includes
finite wide-angle gamma/emitter pairs that vanish in the collinear map.
The older vector photon-pair coefficient uses a different projection,
`y=2*k_gamma.p_selected/M^2`. Its projected scalar kernel therefore differs.
Tests compare the COMPLETE new photon projection plus geometry against that
independently derived coefficient for three pair-mass moments; comparing
only the two projected kernels would be incorrect.

Hard sources use a signed energy-weighted linear-in-z projection, with exact
whole-plus endpoint action. Zero-padded causal FFT convolution transports
them into S/O/G and their first-order coefficients. Tail cells cannot wrap
into observed support. The original hard M determines log(M^2/mu_frag^2),
including off-shell Z parents migrating into the observed mass window.

## Common coupling, overlap and expansion

Let a0=alpha_MSbar(mu0)/(2*pi). The hard coupling a_h is evolved from that
same input to the ORIGINAL M, held there when mu_frag varies. An implicit
scalar two-loop solution avoids repeated matrix evolution in host events;
it is checked against the independent step-doubled RK4 moment solver.

Let H_D and G_D be the projected and geometric hard actions with full
coupled FFs, H_I and G_I their delta-parent limits, and H_c the continuum
bare projected coefficient. The complete matching is

```
raw     = FF_resummed + a_h*(H_D + G_D)
overlap = a0*(C1_FF_cont + H_c - C1_FF_lattice - H_I)
matched = raw + overlap
order1  = Born + a0*(C1_FF_cont + H_c + G_I).
```

The bare geometry cancels between the continuum and lattice overlaps using
identical real samples and local projections. This restores the continuum
NLO term even on a coarse lattice. If H_first and G_first are the coefficients
of a0 in their measured coupled hard actions, the explicit second order is

```
C2 = C2_FF + H_first + G_first
     + beta0*log(M^2/mu0^2)*(H_I + G_I),  beta0=2/3.
```

beta0 is zero in fixed-coupling mode. Photon conversion contributes already
here. Tests verify cubic remainders after subtracting this C2 for W± and Z. W
energy-moment tests also show that factorisation variations contain no
retained alpha² mass logarithm; removing photon parents fails this check.
Products of finite insertions retained by full D convolutions are specified
higher-order terms; unknown constant two-loop boundaries/hard terms are not
promoted to NNLO accuracy. Unlike the separate photonic LL comparison, this
branch does not impose a single finite insertion globally.

## Host interface

`--coupled-hard --hard-samples N` selects the operator and requires coupled
evolution. `--mu-frag` exposes the factorisation-scale variation without
changing the QCD scales or a_h(M). The earlier 26 columns retain their
meaning; seven new rate-only columns give the projected and geometric hard
contributions, their sum, the complete matched rate, explicit a0^2 C2,
matched-minus-continuum-order1, and the paired PHOTOS-minus-matched difference.
All 33 columns enter the raw iteration audit. Inactive-mode columns are zero
and the audit skips their inapplicable algebraic identities.

MCFM real/dipole entries share the original real event's radiation stream
and are summed before squaring. Full Born angular coverage and decay-independent
production scales preserve the spin-tomography identity. New hard columns
are rates; original-parent eta cannot label a radiative hard spectrum.

PHOTOS remains a diagnostic: its fixed alpha(0), pair-free content and
approximate Z radiation differ from this one-family running MSbar reference.
Equal numerical starting alpha is not a finite scheme conversion or an
aligned physical NLL comparison. Common iteration means/covariances retain
the earlier statistical caveats. See the
[benchmark](../benchmarks/coupled_lab_reference.md) for actual completed runs.

At zero alpha the analytic branch reduces exactly to its massless Born
measurement and PHOTOS reduces to its massive mapped Born measurement.
The small difference between those Born maps is retained in the paired
comparison; it is not radiation and is not forced to zero.

Grid, geometric-cutoff and factorisation-scale runs retain the same QCD
production stream. The audit also pairs equal-call raw iteration means,
checking identical Born controls in every iteration, to estimate the
refinement shift's chunk error. This uses ten chunks per contribution; it
does not recover within-event cross-products or certify rare-tail errors.
The conservative separate-run comparisons remain the required pass/fail
controls. A measured grid shift need not be statistically resolved.
