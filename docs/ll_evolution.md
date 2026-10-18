# Fixed-coupling non-singlet LL reference

`LLNonSinglet` evolves a delta boundary in the restricted photonic/non-singlet
theory. It does not include mixing, running α, finite O(α) boundary conditions,
NLO kernels or hard decay. The first-order matrix FF implementation remains a
separate reference; it has not been multiplied into this LL evolution.

For `t = alpha/(2*pi) log(mu_frag²/mu_0²) ≥ 0`, the continuum solution is

```text
D_N(t) = exp(t P0(N))
P0(N) = 3/2 - 2 [psi(N+1)+gamma_E] + 1/[N(N+1)]
P0(z) = [(1+z²)/(1-z)]+
```

The scalar fixed-coupling exponent and LL boundary follow the evolution-operator
construction in [Bertone et al., sections 4 and 5.2.1](https://arxiv.org/html/1911.12040v3#S5.SS2.SSS1).
The LO non-singlet kernel is common to time-like and space-like evolution; this
does not justify importing NLO PDF kernels into FF evolution. The Mellin kernel
above is also checked against the existing numerical plus-distribution action.

## Numerical representation

Set `u=-log(z)`. Multiplicative momentum fractions become additive log losses,
with real-emission measure

```text
nu(u) du = exp(-u) (1+exp(-2u))/(1-exp(-u)) du.
```

Suppress unresolved jumps `0<u<epsilon`. Integrate this positive measure exactly
in each lattice cell, placing its rate at the cell centre `j*h`, `j≥1`. The first
cell is clipped to epsilon. No jump is placed at zero, so epsilon must be at
least `h/2`. A delta at zero log loss carries the negative total virtual rate.
Exponentiating that generator gives the full multiple-emission solution for the
regulated discrete kernel. This is an all-order exponential, not a finite
perturbative expansion or a large-z-only approximation.

The code uses a local radix-two FFT to exponentiate the periodic convolution
generator. It recalibrates phase factors to limit floating-point drift. The
number-conservation Fourier mode is exactly one. No FFT library or other
external dependency is required. The inverse weights are not clipped or
renormalised; negative, imaginary and normalisation residuals are reported and
rejected if too large.

There is an explicit no-resolved-emission atom
`exp(-t * integrated_resolved_rate)`. It represents the numerical unresolved
endpoint. It is not a physical finite delta in the continuum solution at t>0;
its weight changes as epsilon tends to zero. Both endpoints and resolved
weights participate in observable and two-leg convolution.

## Numerical controls and bounds

There are three distinct approximations:

1. Unresolved energy loss below epsilon. Since `u*nu(u) ≤ 2`, its mean omitted
   log loss is bounded by `2*t*epsilon`.
2. Nearest-centre cell placement. Each resolved jump shifts by at most `h/2`,
   giving an additional mean log-loss bound `t*resolved_rate*h/2`.
3. Finite periodic log span U. The code bounds the probability of crossing U
   by Chernoff's inequality with tilt s=0.9, using the **discrete** jump rates:
   `Pr(S≥U) ≤ exp[t sum_j rate_j (exp(s*j*h)-1) - s*U]`.

The rate for jumps larger than the final cell is separately bounded by its
integrated tail times t. The wrap bound concerns the regulated discrete process;
it does not bound the first two errors. The driver refuses excessive wrap
probability. The mean-loss bound controls Lipschitz measurements in u; for
N≥1 moments its factor is at most N−1. Cut observables need their own refinement
study, as their bin boundaries are discontinuities. No single reported bound is
a theory uncertainty or a complete bound on a cut observable's numerical error.

Numerical moments are restricted to real N≥1 to avoid amplifying small-z FFT
roundoff exponentially. Continuum `exact_moment` accepts real N>0. Evolution is
forward only. The non-singlet lepton N=1 moment is conserved; its N=2 moment
decreases as energy is radiated. The scalar model is not a closed matrix in
which the lepton's N=2 moment alone would sum to one.

## Observable action and checks

`apply` evaluates a measurement at every lattice fraction, including z=1.
For identical FFs on both legs, `decay_ll_acceptance` partitions the double sum
by which z is smaller. For fixed smaller z the angular acceptance is constant
in the other leg; the mass window selects an interval of the other log-loss
indices. Prefix sums evaluate this exactly for the lattice weights without a
quadratic tensor loop. A literal double sum validates this reduction for pT,
eta, mass-only cuts and migration from parents above the mass window.

Checks include identity at zero evolution time, conserved number, the endpoint
atom, numerical versus analytic moments, composition by direct discrete
convolution, first/second perturbative coefficients, a fiducial first-order
expansion, and simultaneous grid/cutoff convergence. An independent Python
implementation integrates the jump measure with adaptive Simpson and sums the
Poisson series by direct convolution on a coarse lattice; it checks endpoint
weights and the fiducial result without using the FFT or its analytic cell
integrals.

`leptonff_ll` prints separate baseline, finer-grid and finer-cutoff predictions.
Their residual dependence is part of the numerical study. Pure LL starts from
the delta at mu_0; starting-scale cancellation from the finite boundary is not
claimed. The output's LL first-order result is the appropriate expansion
benchmark, rather than the separate full-O(α)-FF driver.
