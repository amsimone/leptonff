# Conservative resolved-fraction cells

Version 0.19 adds an optional numerical measurement for the coupled W/Z hard
match. It replaces resolved lattice point masses with uniform signed cells,
retaining the exact z=1 atom. It changes neither the evolution kernels nor
the physical flavour/coupling prescription. The default atomic measure is
retained to reproduce the previous reference exactly.

## Why change the measurement

At a finite grid, an atomic fraction measure gives an extra jump whenever a
node crosses an observed cut. These artificial surfaces enter both the
fragmentation action and the locally subtracted hard geometry. They can add
variance to QCD real/counterevent cancellation and obscure a grid comparison.
The cell measure integrates the resolved cut surfaces continuously. The
physical Born endpoint and its cut remain exact; no claim is made that every
source of the hard-estimator variance is removed.

This is a different finite-grid discretisation, not a new physics model and
not an unbiased re-estimation of the old finite-grid point masses. Both
prescriptions must be compared under grid refinement. Their complete
continuum O(alpha) term is identical by overlap restoration.

## Signed measure and conservation

Let z_j=exp(-j*h), d=tanh(h/2), and let w_j be a signed NUMBER weight, obtained
from the energy-tilted evolution only on positive measurement support.
For j>0 use

```
l_j = z_j*(1-d),   u_j = z_j*(1+d)
dN_j(z) = w_j/(u_j-l_j) dz,   l_j<z<u_j.
```

Adjacent resolved cells meet at their harmonic-mean edge. Each cell is
centred at its original z_j, so its integrated number is w_j and energy is
w_j*z_j. The endpoint is w_0*delta(1-z), without broadening. Negative weights
are retained. There is no clipping, fitted correction or number
renormalisation. Global singlet number integrals need not exist.

One extra node is retained when its cell intersects the positive support
boundary. Queries below minimum_z are rejected, including in inverse-fraction
prefixes. The full hard/evolution small-z tail is never folded into retained
cells. The evolution's weighted wrap/tail controls remain separate from the
measurement's discretisation error. A cell intersecting a lower cut q has
its node at z_j>=q/(1+d); the number-measure test on energy weights is bounded
by (1+d)/q. For a two-leg measurement the corresponding bounds multiply.

Higher moments change: a resolved cell's second z moment is
`w_j*z_j^2*(1+d^2/3)`. Preserving number and energy does not make all moments
exact. An independent smooth D(z)=1 fixture shows second-order CDF convergence;
that does not establish a second-order error bound for the physical endpoint
singularities or a certified hadronic error.

## Exact two-leg cuts

Define signed tails

```
C(q) = integral_(z>q) dN(z)
V(q) = integral_(z>q) dN(z)/z.
```

They are piecewise linear and logarithmic on resolved cells, with the exact
endpoint atom. Prefix evaluation uses one binary search for both tails.
For a first-leg cell segment a<z1<b and a second-leg pT lower bound l2, put
L=max(l2,K/b), U=max(l2,K/a). Reversing the integration order gives

```
T_K(a,b;l2) = integral_a^b C(max(l2,K/z1)) dz1
            = (b-a)*C(U) + b*(C(L)-C(U)) - K*(V(L)-V(U)).
```

The exact rectangle action for `mass_lo<z1*z2<mass_hi` is
`density1*(T_mass_lo-T_mass_hi)`. The first-leg atom is evaluated separately
with a strict interval query. The formula also handles a second-leg atom and
mixed atomic/cell measures. Zero-weight cells, empty mass support and exact
atom-only partners have fast paths. No inner numerical quadrature or random
fragmentation fractions enter the host callback.

## Matching and host interface

`SignedFractionTable`, coupled `LabFragmentation`, and `CoupledLabDressing`
accept `FractionMeasure::conservative_cells`. The same choice acts on:

- S/O/G and their measured real and collinear hard events.
- C1/C2 and the Born-parent FF overlap.
- The bare projected hard source, its transported sources and their C1.

The [coupled match](coupled_lab_matching.md) therefore retains its continuum
NLO restoration and explicit second-order coefficient. W±/Z expansion tests
pass for both measures, leaving cubic remainders after subtracting C2. Using
cells in only the resummed term, while subtracting the old atomic overlap,
would spoil this property and is not supported by the provided runner.

`--fraction-cells` requires `--coupled-hard --hard-samples N`. All 33 rate
columns and the original-real radiation streams are retained. The flag and
an active-callback log marker are recorded in every manifest. The fixed-alpha
photonic branch is unchanged.

## Verification and interpretation

Independent tests integrate signed cell overlaps, reciprocal moments and
hyperbolic rectangle areas directly, including negative weights and both
strict endpoint atoms. They check number/energy conservation, resolved CDF
continuity, mixed-measure symmetry, positive support, and convergence to a
separately specified smooth continuum density. The coupled tests check the
actual cell-measured C1/C2 expansion and zero coupling.

Host comparisons use the same QCD streams for cells and atoms, common raw
iteration means and directly accumulated higher-order differences. Observed
variance ratios are sample diagnostics; they are not guarantees over random
seeds, phase space or observables. CPU-time ratios are reported separately.
Independent-seed checks, cutoff/subtraction variations and grid shifts are
retained in the [benchmark](../benchmarks/fraction_cell_reference.md).
Physical flavours/input conversion, hard differential spectra and a complete
physics/numerical/statistical uncertainty budget remain outstanding.
