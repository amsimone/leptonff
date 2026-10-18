# Coupled fiducial hard matching

Version 0.11 joins the complete one-family NLO time-like FF evolution with all
O(alpha) hard-parent channels for the fiducial unpolarised vector decay at rest.
`CoupledFiducialDecay` is the reusable hard operator; `leptonff_coupled` provides
the reference rate, numerical refinements and fragmentation-scale scan.

This is an O(alpha) hard-matched NLL reference in the **closed one-family MSbar
model**, with symmetric pT, eta and opposite-charge pair-mass cuts. It supplies
an inclusive identified-pair density, not an exclusive two-muon event rate.
Physical Z spin/axial terms, boosts, active flavours, thresholds, hadronic
vacuum polarisation and a production host are not supplied by this increment.

## Measurement and channels

Use the [pair-counting contract](photon_pair_matching.md): one observed lepton
per selected hard direction, with no veto on extra pairs. Let S, O and G denote
`D[mu- -> mu-]`, `D[mu- -> mu+]` and `D[gamma -> mu-]`, respectively. Charge
conjugation supplies the other entries. The three parent-pair channels are
`mu- mu+`, `gamma mu-`, and `gamma mu+`. With B the Born cut functional,

```text
raw = B[S,S] + B[O,O]
    + a_h * (Hm[S,S]+Hm[O,O]+2*Hg[G,S+O]),
a_h = alpha_MSbar(M)/(2*pi).
```

Hm includes its virtual/plus mass projection and finite real geometry; Hg
includes the subtracted photon-primary mass projection and finite geometry.
The factor two sums primary charges, not a second duplicate sum over observed
charge assignments. The photon term begins at alpha²*log(M²/m²) and is needed
at NLL. The [moment reference](photon_pair_matching.md) validates the same
channel contraction and common hard coupling without fiducial cuts.

Both FF running and a_h use the same explicit MSbar input at mu_0. During a
mu_frag variation the hard coupling remains at M. The reference config's alpha
is interpreted as MSbar alpha(mu_0), without an on-shell conversion. Two-loop
running and the complete NLO matrix are used. Products/exponents keep selected
terms beyond NLL; no global single-boundary restriction from the older scalar
LL driver applies. Constant two-loop FF boundaries and H2 are not supplied.

## Hard mass projections

The source FF grid stores energy-weighted measures `z*D(z) dz`, with
`z_j=exp(-j*h)`. A pair passing the cuts necessarily has each fragmentation
fraction at least

```text
minimum = max(2*pt_min/M, (mass_min/M)²).
```

The operator requires mass_min>0, nonempty support and a grid extending below
minimum. It retains only the causal prefix above minimum. The hard distributions
are projected with the same linear-in-z basis as the FF generator. In the first
cell, with `z1=exp(-h)`, endpoint and neighbour weights are

```text
w0 = H[ z*(z-z1)/(1-z1) * theta(z-z1) ],
w1_first_cell = H[ z*(1-z)/(1-z1) * theta(z-z1) ].
```

The full distribution action in w0 includes the entire virtual/plus subtraction
below z1, including outside the retained prefix. Every other retained node
receives both adjacent cells; its lower adjacent cell is included even when
the next node lies outside support. No unrepresented real tail is moved into
the endpoint, and there is no imposed normalisation.

Use ordinary causal convolution `C_j=sum_(k=0)^j H_k*D_(j-k)`. Undo the energy
tilt only on the retained support before the Born bilinear measurement. This
hard convolution has neither periodic wrap nor a separate hard tail regulator.
The underlying FFT FF evolution still has its reported weighted wrap/tail/step
diagnostics. Those diagnostics are not a bound on total observable error.

For reuse across scale scans, project Hm/Hg at mu_frag=M and their LO kernels.
The scale terms are exactly `Hm(mu_frag)=Hm(M)+2*log(M²/mu_frag²)*P_mu_mu` and
`Hg(mu_frag)=Hg(M)+log(M²/mu_frag²)*P_mu_gamma`.

## Finite geometry and signed FF quadrature

The muon geometry is the already validated locally finite real-minus-collinear
projection used by the scalar LL driver, now applied to each required coupled
FF pair. The photon geometry is the fixed-fraction coefficient validated in
[the photon-parent increment](photon_pair_matching.md).

A dense FF tensor inside each two-dimensional geometry integral would be costly.
Instead, panel the resolved support uniformly in sqrt(1-z). Each panel has four
fixed Gauss abscissae in z. Deposit each resolved lattice weight using the cubic
Lagrange basis; retain the endpoint atom separately. Empty panels are removed.
The weights are signed, with no clipping. Gauss abscissae only choose interpolation
locations here: this is not a positive Gaussian quadrature of a probability.

Crucially, this projection is **linear in its input weights**. The raw FFs,
first-order coefficients and second-order coefficients therefore use the same
measurement operator. A density-dependent compression would not automatically
have the expansion used to subtract its overlap. Number and z powers through
three of each truncated unweighted measure are preserved to roundoff; those
identities do not bound interpolation errors for sharp cuts. Compare panel
counts and hard quadrature/sliver settings for an actual precision statement.

Geometry matrices are precomputed once, independent of the FF coupling/scale.
Muon geometry is symmetric in its two fractions; photon geometry is ordered.
The numerical slivers omit `1-z<soft_cutoff` and `1-|c|<collinear_cutoff` in the
muon remainder, and `t/(1-y)<collinear_cutoff` in the photon remainder. They
regulate already finite integrals, not the MSbar subtraction. Positive pair mass
excludes soft photons in the photon-primary channel exactly.

## Complete overlap matching and expansion

The matcher extracts the **actual numerical** C1 of the combined FF and hard
operator, then restores the continuum full O(alpha) coefficient:

```text
matched = raw + a0*(C1_continuum-C1_lattice),
a0 = alpha_MSbar(mu_0)/(2*pi).
```

C1_lattice includes the two FF legs and the hard endpoint/mass projection and
finite geometry. Its subtraction is not limited to a collinear logarithm.
C1_continuum is the existing fiducial massless hard+FF result, with identical
hard geometry controls. This guarantees the established leading-power O(alpha)
target subject to the reported integration/sliver errors.

The returned coefficient2 is the explicit alpha(mu_0)² coefficient of this
numerical prescription. It contains both second-order FF legs, products of
first-order FFs, Hm acting on either first-order leg, the hard-running term
`beta0*log(M²/mu_0²)*Hm[I,I]`, and the photon term `2*Hg[C1_gamma,I]`.
The overlap restoration is linear in a0 and leaves coefficient2 unchanged.

Tests show a cubic remainder after removing Born, a0*C1 and a0²*C2, and verify
the common hard coupling against the independent moment interface. Under a
mu_frag variation C1 is invariant and the remaining C2 scale shift has no
retained mass logarithm. Constant NNLL scale dependence remains. Grid and
quadrature artifacts in C2 must be resolved numerically; missing H2/d2 terms
are not supplied by additive matching.

Python independently builds source cells in z, evolves a causal power series
without FFT/RK steps, reconstructs hard endpoints from analytic moments, and
checks literal Born/hard contractions. Geometry is integrated with allowed
azimuth intervals instead of spherical-cap intersections. These are independent
algorithm/current checks within the same model, not external end-to-end physical
NLL validation. See the [benchmark](../benchmarks/coupled_decay_reference.md).

## Reproduction

```sh
./build/leptonff_coupled configs/reference.conf
./build/leptonff_coupled configs/reference.conf 8192 8 --scan \
  > results/coupled_decay_fine_scan.json
python3 python/check_coupled_decay_reference.py build/leptonff_coupled \
  configs/reference.conf --input results/coupled_decay_fine_scan.json \
  --output results/coupled_decay_reference.json
```

`--dump` exports the causal hard prefixes, retained FF/coefficient arrays,
fixed fractions and geometry matrices for a single point. `--scan` doubles
FF panels, then grid resolution, then span at fixed finest spacing; next it
divides both hard slivers and geometry quadrature tolerances by four. The last
two points vary mu_frag/M by 1/2 and 2 with the same final operator. Top-level
hard controls refer to the baseline; the `hard_slivers` and scale points use
the controls divided by four. Coarse independent checks use 256 FF bins/two
panels and the unrefined controls. Saved manifests identify code, executable,
configuration and any attached production scan.
