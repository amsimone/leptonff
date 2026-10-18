# Evolved O(α) boundary and FF-only matching

Version 0.3 adds the signed MS̄ diagonal boundary to the non-singlet LL operator.
This is a step toward NLL bookkeeping, not a complete NLL or hard-decay result.
The boundary is the existing coefficient `d1(mu_0)` documented in
[kernel_sources.md](kernel_sources.md); it is not an extra photon shower.

## Order convention

Define `a=alpha/(2*pi)`, `L=log(mu_frag²/mu_0²)`, `U=exp(a L P0)` and
`B=U⊗d1`. A single evolved FF is represented as `U+aB`. For two muons use

```text
R = (U ⊗ U)[F] + a (B ⊗ U + U ⊗ B)[F].
```

There is at most one boundary insertion across both legs. Multiplying two
unrestricted `U+aB` arrays would also add `a² B⊗B`, which is deliberately excluded.
The boundary is neither exponentiated nor constrained to be positive.

The single-leg expansion contains

```text
I + a(d1 + L P0) + a²[L P0⊗d1 + L² P0⊗P0/2] + ...
```

The missing `a² L P1` term is required for NLL; including the boundary-generated
term alone is insufficient. The implementation also omits running, mixing,
physical pair prescriptions and all O(α) hard coefficients.

## Lattice representation and signed endpoint

Use the same cutoff, log-loss cells and phase-space convention as the LL solver.
The boundary's resolved measure is

```text
nu(u) [log(mu_0²/m_mu²) - 1 - 2 log(1-exp(-u))] du.
```

The cell integrals are analytic. Their calculation uses stable differences
after substituting `s=1-exp(-u)`; it avoids subtracting nearly equal primitive
values at large u. The zero-loss coefficient is minus the sum of resolved
coefficients, implementing the whole-plus prescription. The number eigenvalue
is zero, not one. The FFT convolves this signed boundary with U linearly.

Initial/evolved number residuals and imaginary roundoff are reported; real
negative coefficients remain part of the result. A wrap bound for this signed
convolution uses the LL Chernoff bound multiplied by the tilted total variation
of the initial boundary. The omitted large-u boundary tail is separately bounded.
These bounds concern the finite-domain treatment, not soft-cutoff/grid errors.

The two-leg prefix reduction is bilinear in two independent signed arrays. It
therefore handles U×B directly, including negative endpoint coefficients and
the common invariant-mass cut. Tests compare this with a literal tensor sum.

## Complete first-order overlap at finite resolution

The regulated lattice has its own LO generator K and boundary b. Its exact
O(α) coefficient for the decay measurement is

```text
C_lattice = ((b+L K) ⊗ I + I ⊗ (b+L K))[F].
```

The continuum `C_FF` comes from the independently implemented full O(α) FF
coefficient, including its finite terms. The additive result is

```text
R_matched = R + a (C_FF - C_lattice).
```

This removes the **actual** lattice overlap, rather than a guessed logarithmic
term or the continuum LL coefficient. It reproduces the complete FF-only O(α)
result already at finite grid/cutoff. Higher orders still carry numerical
regulator effects and require refinement.

This is not physical NLO matching: the hard O(α) decay is absent. Once a validated
hard decay is implemented, it must enter both the resummed approximation and
its corresponding fixed-order overlap. A physical massive benchmark remains
necessary; it cannot be replaced by agreement of the two FF implementations.

## Validation

- Analytic cell integrals versus existing plus-distribution action on cell cuts.
- Signed FFT convolution versus direct convolution, and number conservation.
- Evolved moments versus `U_N*d1_N`, with regulator refinement.
- Matched-minus-fixed-order remainder scaling as α² under a coupling scan.
- Equal evolution scales giving exactly the continuum O(α) FF result.
- Explicit exclusion of double-boundary terms across two legs.
- Independent Python boundary integration with Simpson, direct signed
  convolution, independent fixed-order coefficient and overlap calculation.

`leptonff_match` uses the same configuration and optional grid arguments as
`leptonff_ll`. It reports every matching component at baseline, doubled grid
resolution and halved cutoff on the finer grid. The refinement shifts are
numerical diagnostics rather than perturbative uncertainty estimates.
