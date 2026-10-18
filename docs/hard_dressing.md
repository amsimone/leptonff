# LL dressing of the finite hard decay

Version 0.6 supplies the missing H1 convolution for the unpolarised vector
parent at rest. Both muons of the real hard event receive independent pure-LL
fragmentation. Together with the existing single evolved boundary, this
generates hard- and boundary-induced subleading logarithmic terms. It still
does not contain the NLO splitting kernel, running, mixing, DY spin density,
boosted parent or physical axial mass effects.

## Perturbative contract

The evolved-FF/subtracted-hard framework is discussed in
[Frixione, 1909.03886v2, section 7](https://arxiv.org/html/1909.03886v2).
The single-finite-insertion formula below is this implementation's restricted
expansion of that framework.

With U the fixed-alpha non-singlet LL solution and a=alpha/(2*pi), retain

```text
H0[U x U] + a*H0[(U*d1) x U + U x (U*d1)] + a*H1[U x U].
```

There is at most one finite insertion globally. In particular, H1 is dressed
with U alone: no H1*d1 or H1*H1 product is introduced. The previous additive
driver remains available as a comparison. Setting the evolution time to zero
gives the complete O(alpha) result, including its finite boundary and hard term.

## Collinear projection

The [fiducial hard coefficient](fiducial_vector_matching.md) splits into the
mass-spectrum distribution applied to a collinear measurement and a finite
geometry remainder. Integrate the mass-spectrum kernel over each LL loss cell;
store the signed endpoint coefficient as 3/2 minus the sum of resolved cells.
Analytic primitives use stable small-z and soft expressions, including an
explicit subtraction of the primitive at z=0 to avoid cancellation in far-tail
cells. Independent numerical integration checks the cell coefficients.

One signed FFT convolution evolves this hard distribution on the first leg.
The existing bilinear prefix-sum measurement supplies U on the second leg.
Symmetry of the identical muon FFs accounts for both collinear endpoints; no
extra factor of two multiplies the already two-leg mass-spectrum coefficient.
An angular-inclusive mass window reduces independently to H1*U(2*t).

The implementation reports number and imaginary FFT residuals. A tilted
variation of the signed hard weights multiplies the two-leg LL Chernoff bound;
an analytic tail bound covers the omitted hard measure beyond the grid.
The rate-level wrap/tail bound is checked against the grid tolerance.

## Finite geometry with two fragmented muons

For real kinematics (z,c), replace p1,p2 by f1*p1,f2*p2. The directions stay
fixed and the observed pair mass becomes M*sqrt(z*f1*f2). Recompute both pT
bounds before averaging over the common rotation. For c>0 the local collinear
projection is B12=isotropic_acceptance(M,z*f1,f2). The reflected angular
region uses B21 with the two fractions exchanged. Both terms are required:
multiplying an undressed correction by a scalar acceptance would miss this
geometry and would not implement the hard convolution.

The finite remainder uses the same soft and angular slivers as version 0.5.
For each pair of FF nodes the integration follows the observed mass bounds and
splits at the pT/eta transitions. Swapping the node fractions handles the
reflected c<0 region. This avoids integrating a singular real contribution and
its collinear counterterm separately.

The full LL lattice has too many entries for a literal tensor sum at every
hard phase-space point. Retain its endpoint atom exactly and group the positive
resolved measure in sqrt(1-z) panels. A two-node Gaussian rule per panel
reproduces its first four polynomial moments (powers 0 through 3); centred
moments keep the Jacobi construction stable near z=1. Signed FFT roundoff is
retained in a separate weighted node, preserving its number and first z moment.
The implementation checks the resulting moments against the uncompressed
lattice. This quadrature is a numerical approximation: node-count refinement
is necessary, especially for sharp cuts. It is distinct from LL grid and
soft-cutoff refinement.

## Matching and checks

Let P_lat be the projected hard lattice measurement and G_Q the finite
geometry evaluated with the reduced LL quadrature. The coefficient used is

```text
H1_dressed = H1_continuum
             + (P_lat[U,U] - P_lat[I,I])
             + (G_Q[U,U] - G[I,I]).
matched = FF_matched + a*H1_dressed.
```

The zero-evolution hard term is subtracted with the same numerical controls
before restoring the continuum fixed-order coefficient. This prevents lattice,
sliver or quadrature errors from spoiling O(alpha) matching. The old FF matcher
already removes its entire discrete O(alpha) overlap. Cached hard and evolution
configurations must agree.

Tests compare signed FFTs to literal cyclic sums, the hard-cell primitive to
distribution actions, compressed moments to the original lattice, and both
local collinear measurements to their analytic limits. Mass-only tests use the
independent U(2*t) reduction. Coupling scans establish the alpha-squared onset
of both the complete matched remainder and the change from additive matching.
No-evolution and inclusive-number tests check that finite products have not
been introduced.

Python independently integrates jump and hard cell rates, sums a Poisson
series without FFTs, and directly convolves signed arrays. Its angular check
uses direct azimuthal quadrature rather than the spherical-cap formula. The
independent angular fixture is deliberately coarse enough for a full check;
the production calculation separately reports quadrature, grid, LL cutoff and
hard-geometry refinements. See [the benchmark](../benchmarks/hard_dressed_reference.md).

The H1 dressing gap is closed for this restricted model. Full NLL remains
gated on the time-like P1 source/mapping audit and the coupling/flavour content.
The physical spin-correlated Z and host-adapter gates remain open.
