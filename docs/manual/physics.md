# Scientific contract

## Fragmentation conventions

The core uses a **parent to observed** index convention. The one-family basis
is `(muon_minus, muon_plus, photon)`. Kernels are timelike; PDF/spacelike kernels
are not interchangeable with them. The fixed-order boundary and kernels use
MSbar factorization and the expansion parameter $a=\alpha/(2\pi)$.

For one family, the first-order fragmentation function is

$$
D_{i\to j}(z,\mu_F)=\delta_{ij}\delta(1-z)
+a\left[d^{(1)}_{i\to j}(z,\mu_0)
+\log\frac{\mu_F^2}{\mu_0^2}P^{(0)}_{i\to j}(z)\right].
$$

`fragmentation_coefficient_order1` returns the coefficient of $a$;
`fragmentation_order1` returns the full expression through first order.
The finite first-order boundary is retained.

`Distribution` represents a delta contribution, a regular part and whole-plus
distributions. In particular,

$$
\int_0^1 dz\,[g(z)]_+ f(z)=\int_0^1 dz\,g(z)[f(z)-f(1)].
$$

The prescription acts on the complete numerator and denominator in `g`.
`resolved_density(z)` supplies only the interior real-emission density.
Moments, bin integrals and convolutions include endpoint subtraction terms.

`convolve(f,x)` uses the standard spectral convolution with measure $dz/z$.
`apply(measurement)` instead integrates a fragmentation function against an
event measurement. When a host supplies weighted hard events and the
measurement evaluates rescaled momenta, no additional $1/z$ belongs in that
event action.

## Evolution and accuracy labels

The coupled backend solves

$$
\frac{dE}{d\log\mu^2}=(aP^{(0)}+a^2P^{(1)})E,
\qquad D=E\,[I+a(\mu_0)d^{(1)}].
$$

Kernel order, logarithmic evolution accuracy, boundary order, hard matching
order and production QCD order are distinct settings. A coupled NLL
fragmentation result alone does not establish NLL accuracy for an arbitrary
fiducial observable. Unknown constant two-loop boundary terms are omitted.

`LLNonSinglet` supplies restricted fixed-alpha non-singlet LL evolution.
`LLBoundaryNonSinglet` adds a single finite boundary insertion. The conditioned
photonic NLL implementation belongs to the matched decay prescription below.

## Coupling and flavour inputs

`EvolutionTheory` separates physical settings from solver/grid controls.
`AlphaInput::msbar(value, scale)` requires that scale to equal `initial_scale`.
The interface rejects a mismatch; callers must run the coupling explicitly.
`AlphaInput::on_shell(value)` performs the implemented one-loop leptonic
conversion at the initial scale. Hadronic conversion requires separate sourced
input and an explicit matching prescription.

`QEDCouplingMode` independently selects fixed, one-loop or two-loop running.
The optional `QEDFlavourScheme::leptons()` activates electron/muon/tau thresholds
and the implemented finite matching. An empty flavour scheme selects the
legacy one-family backend. With additional active flavours, the returned
3-by-3 muon/photon block is a subblock; momentum conservation concerns the full
flavour system. Extra identified leptons describe inclusive densities, without
an exclusive extra-pair veto.

Advanced pointlike-quark/HVP and LO QCD study settings remain in
`MomentEvolutionParameters`. They do not constitute fitted hadronic fragmentation
functions or complete mixed QCD/QED evolution.

## Moment and x-space domains

Mellin moments use $D(N)=\int_0^1 dz\,z^{N-1}D(z)$. The coupled real-moment
solver accepts $N>1$. Individual channels can have divergent number moments;
the one-family fixed-order interface checks each channel's domain.

`QEDXSpaceEvolution` stores **$zD(z)\,dz$**, using an energy-weighted log grid.
Its lattice moments require $N\geq2$. Actions undo the energy tilt only on
positive lower support; sharp bins need grid convergence checks. Small-z tails,
periodic wrap and integration/step errors are distinct diagnostics. Signed
weights are retained.

## Matched W/Z prescription

The matched module combines fragmentation with finite spin-correlated W/Z
decay responses. Photonic NLL uses the explicitly conditioned shape
$E_{\rm raw}(N)/E_{\rm raw}(1)$, consistent with a supplied inclusive branching
normalization. The published source kernel's nonzero number moment remains
recorded in diagnostics. This shape is distinct from a raw MSbar fragmentation
function or an absolute NNLO partial width.

Across multiple decay legs, NLL dressing retains **one subleading insertion in
total**: boundary, hard coefficient, or the NLO photonic evolution insertion.
Multiplying independently truncated leg factors can introduce unwanted terms.

The Z hard response supports parity-even, charge-symmetric measurements.
Complete off-shell EW corrections, complete mixed QCD/EW corrections and native
DYTurbo NNLO+NNLL predictions are outside the implemented scope.

## Production host contract

A host supplies compatible pre-FSR production points or spin density tensors.
Tensor axes and recoil conventions must match; the Drell–Yan helper explicitly
rotates Collins–Soper angular weights into the library's Cartesian convention.
Projected hard matching needs full Born angular coverage, without original
lepton-dependent cuts or scales that invalidate that projection.

Keep signed real/subtraction groups intact, use a common radiation seed for
their mapped points, and accumulate group covariance. A production archive
stores weighted points; the current MCFM exporter supplies central QCD weights.
It is not a PDF-independent interpolation grid. Already integrated histograms
do not generally retain the information needed to reconstruct FSR migration.
