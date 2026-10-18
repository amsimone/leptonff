# Coupled x-space evolution and endpoint control

Version 0.9 supplies `QEDXSpaceEvolution`, a signed log-loss lattice for the
complete one-family MS̄ time-like NLO matrix. It uses the same kernels, running
coupling, initial boundary and matrix ordering as the
[moment reference](moment_evolution.md). The two-loop input is explicitly
`alpha_msbar_at_mu0`. No coupling conversion or flavour threshold is implicit.

## Energy-weighted measure

The stored arrays represent `z D(z) dz`, rather than `D(z) dz`. This energy
tilt makes the singlet and photon small-z terms integrable. Lattice nodes are
`u_j=j*h`, `z_j=exp(-u_j)`, `h=log_span/bins`, including the endpoint j=0.
Multiplication of fragmentation fractions becomes addition of losses, and
the energy weights multiply consistently under convolution.

The public `weights`, `propagator`, `kernel0`, `kernel1`, `boundary1` and
expansion-coefficient arrays all use this same tilt. `moment(N)` sums
`weights[j]*z_j^(N-2)` and requires N>=2. A measurement action requires an
explicit positive lower bound; it undoes the tilt only on that support.
Full inclusive identified-lepton number integrals need not exist, and dividing
every far-tail FFT cell by its tiny z would amplify floating-point noise.

## Endpoint projection

Project a test function linearly in **z**, between consecutive nodes. For
each interior segment, integrate `z*P(z)` times the two interpolation weights
and deposit the results at the segment's bounding nodes. In the first segment
the neighbour test is `z*(1-z)/(1-z_1)`, which vanishes at z=1. Applying the
full distribution to this test gives a finite neighbour coefficient including
the regular and plus terms. The endpoint coefficient is then

```text
K[0] = P(N=2) - sum_(j>0) K[j] - integral_0^z_min dz z P_real(z)
z_min = exp(-(bins-1)*h)
```

This identity retains the source delta and all virtual subtractions. It uses
each individual distribution's N=2 moment; no parent row is renormalised.
The omitted real tail is subtracted explicitly, so it is not returned to z=1.
The same construction applies to the signed O(alpha) boundary.

The N=2 kernel moment is fixed by construction. N=3 is also reproduced because
the additional test z is linear in the interpolation coordinate. N=5,8 and
the fiducial measurements provide independent projection-convergence tests.
There is no separate soft cutoff: h controls the first-cell resolution. Sharp
cut errors need not converge at the smooth-test second-order rate.

## Coupled evolution, bounds and diagnostics

FFT each weighted kernel and boundary. For every frequency solve
`dE/dL=(a*K0+a²*K1)E` with the existing step-doubled matrix RK4 method, then
form `D=E*(I+a_initial*B1)`. Hermitian frequency symmetry gives real signed
weights after the inverse FFT. Neither Fourier zero modes nor weights are
adjusted to enforce conservation. Negative weight mass and imaginary residuals
are reported separately.

To control cyclic aliasing, construct a positive comparison generator. Retain
each signed diagonal zero-loss entry, take absolute off-diagonal zero-loss
entries, and sum absolute resolved cells weighted by `exp(sigma*u_j)`.
Its maximum row sum is `mu_r(sigma)` for kernel r. Let `B(sigma)` be the
absolute boundary row norm, including its identity term. With `I1=integral a`,
`I2=integral a²` and `0<sigma<1`, the weighted wrap comparison is

```text
wrap <= B(sigma) exp(I1*mu_0(sigma) + I2*mu_1(sigma) - sigma*log_span)
```

The restriction sigma<1 matches the small-z domain of the energy-weighted
matrix. Absolute tail integrals `T0,T1,TB` give a separate perturbation bound:

```text
R = I1*T0 + I2*T1
tail <= exp(I1*mu_0(0)+I2*mu_1(0))
        * (B(0)*expm1(R) + a_initial*TB*exp(R))
```

The tail quadrature's estimated error is added to its absolute integral.
These comparison bounds concern weighted wrap/tail effects; they do not
bound cell projection, quadrature or floating-point errors in the full result.
Too-short domains are rejected. For a bounded measurement supported on
z>=z_cut, undoing the tilt multiplies a weighted variation bound by at most
`sup(abs(F))/z_cut`.

The frequency step-refinement diagnostic also reports the conservative
`3*sqrt(bins)*max_frequency_shift*B(0)` row-variation bound on the compared
step grids, using Parseval and the boundary convolution norm. This is a
refinement diagnostic, not a proof of the remaining integration error.

## FF fiducial validation and remaining hard matching

`decay_xspace_ff_match` applies the Born two-direction measurement to
`D_same*D_same + D_opposite*D_opposite`: one identified opposite-charge lepton
per direction, counting inclusive multiplicities. It is not an exactly-two-muon
event selection. Undo the tilt only above the physical pT/mass support.
Subtract the actual lattice O(alpha) overlap and restore the continuum FF
coefficient, preserving the complete first-order FF result. The two legs retain
their boundary products, as in the NLO-generator moment prescription; this
includes some terms beyond NLL, unlike the older globally linear LL reference.

This calculation has **no hard H1 contribution** and is not a full fiducial
NLL cross section. Connecting the full-QED evolution to the matched decay
requires the hard photon-parent contributions at the retained order, consistent
hard-coupling conversion and an explicit pair-selection prescription. Adding
the old muon-only hard coefficient would not complete those requirements.

The independent Python check integrates cells in z using fixed Gauss rules
and evolves them with a causal, non-periodic power series in the evolution
interval. It uses neither FFTs nor RK4. It checks every channel/cell, first and
second expansion coefficients and a literal fiducial pair sum. See the
[benchmark](../benchmarks/xspace_reference.md) for the achieved grid/domain
refinements and reproduction commands.
