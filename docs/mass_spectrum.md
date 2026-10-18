# Matched bare-pair mass spectrum

Version 0.12 adds `CoupledMassSpectrum` and `leptonff_mass_spectrum`. They return
mass-bin integrals of the [coupled fiducial match](coupled_fiducial_matching.md),
including every supplied hard-parent channel and the complete O(alpha) overlap.
Each point contains its fixed-order target, raw and matched result, FF/hard
components, first/second expansion coefficients, total and sum over bins.

The observable is the inclusive opposite-charge **pair mass**, with the original
symmetric pT/eta cuts. It is not a single-muon pT/eta spectrum. The model remains
a pure unpolarised vector at rest, one lepton family, explicit MSbar input,
inclusive distinct-parent pairs and no extra-pair veto. The parent has a fixed
mass M: bins above M are exactly empty, rather than a physical off-shell DY tail.

## Shared cumulative measurement

Edges must be finite, strictly increasing, and exactly partition the outer
fiducial mass window. Keep the entire outer cut contract, including positive
mass_min and rejection of a **strict outer upper cut exactly at M**.

For each edge below M construct the cumulative rate with pair mass above that
edge, retaining the same outer upper cut. All cumulatives share the original
FF support, hard source projections and cubic interpolation fractions. Only
the measurement and its finite-geometry matrices change. Cumulatives at or
above M, and at the final outer edge, are identically zero. A bin is the
difference between its two adjacent cumulatives.

This is a linear partition of each additive term, not a normalisation imposed
after computing independent bin operators. Recomputing a different FF
compression for each bin would instead change the numerical measurement.
`restricted_mass_window` exposes the shared-support restriction for other audits.

Internal bins have convention `(lower,upper]`, with the original strict outer
cuts retained. If an internal edge is exactly M, the bin below it owns the Born
atom and its arbitrarily soft neighbourhood; the bin above it is zero. This is
infrared safe. It differs from excluding Born with a strict upper cut at M while
accepting radiation just below M. Floating-point grid-edge ownership is handled
by the common cumulative rule, so no bin is silently dropped between two strict
independent cuts.

Couplings are shared metadata, and compression diagnostics are common operator
diagnostics. They are not differenced or summed. All additive rate components
and C1/C2 coefficients are partitioned. The returned total is computed by the
original outer-window operator; it is not obtained by renormalising the bins.

## Conservative mass-cell integration

An integrated fiducial rate can hide bin-edge aliasing. On a log-loss grid,
`f_i=exp(-i*h)`, the Born and projected hard measurements have pair mass
`M*exp(-(i+j)*h/2)`. Point assignment moves a whole signed shell between adjacent
bins as the grid changes, even when the total changes very little.

The spectrum therefore integrates each **resolved total-loss shell** centred at
`u=(i+j)*h>0` over `[u-h/2,u+h/2]`, with constant weight per unit u. For a mass
bin (lo,hi], its fraction is

```text
length( [u-h/2,u+h/2] intersect [2*log(M/hi),2*log(M/lo)] ) / h.
```

The zero-loss shell stays an explicit atom: it is never spread across bins.
Angular cuts and energy-tilt factors are evaluated at the grid fractions as
before. The efficient production implementation uses prefix sums and at most
two transition nodes; tests compare a distinct literal signed tensor/cell loop.

This is a numerical integration prescription whose cell width vanishes under
refinement, not detector resolution or a physical photon/muon recombination.
It does not alter the FF evolution or its moment definitions. The spectrum
selects `CoupledDecayOptions.average_mass_cells=true`; the older rate driver
retains point assignment by default. Both converge to the same observable;
finite-grid totals need not coincide exactly. The source support includes an
extra half-cell in total loss when the mass cut controls the truncation.

The cell rule is applied consistently to the raw FF term, causal hard
projections and the numerical first/second expansion. Finite geometry remains
continuously integrated at fixed FF quadrature fractions. Restoring the full
continuum O(alpha) overlap gives the same fixed-order bin target for either
discretisation. Tests verify that the residual alpha² scale shift has no
retained mass logarithm in each cell-integrated bin.

Closure is necessary but does not establish numerical convergence. In particular,
opposite shifts in neighbouring endpoint bins can cancel in the total. The
driver reports grid/panel/sliver variations per bin; never infer their errors
from the total-rate refinement alone. Bins much narrower than the grid's mass
spacing require their own refinement study.

## Correlated variations and units

Every scale point applies one common variation to the entire spectrum. During
mu_frag/M={0.5,1,2}, FFs and hard logarithms vary together, while alpha_MSbar(M)
stays fixed. Starting-scale variations use mu_0 factors {0.5,2}, transport the
MSbar coupling from the nominal mu_0 to the new value, and recompute the FF
boundaries. The same hard geometry operator is reused. All bins still sum to
the varied fiducial rate.

Starting-scale variation is defined here in the closed one-family model, with
no physical thresholds. Unknown constant two-loop boundaries/H2 are not added;
these variations probe retained versus omitted terms, not a complete physical
uncertainty. Report numerical shifts, mu_frag and mu_0 variations separately.
The bin-wise extrema of a scale envelope are not one common varied spectrum.
There is no sampling noise or statistical covariance in this deterministic
reference, and no covariance matrix is fabricated.

`integral.matched` is dimensionless relative to the massless inclusive vector
Born normalisation. `matched_density_per_GeV` divides that integral by its bin
width. The photon entries already sum both primary charges. Signed terms are
retained without clipping. Inclusive pair density is not an exclusive decay
probability or an absolute cross section in pb.

## Commands and checks

```sh
./build/leptonff_mass_spectrum configs/reference.conf configs/mass_bins.conf
./build/leptonff_mass_spectrum configs/reference.conf configs/mass_bins.conf \
  8192 4 --scan > results/mass_spectrum_fine_scan.json
python3 python/check_mass_spectrum_reference.py build/leptonff_mass_spectrum \
  configs/reference.conf configs/mass_bins.conf \
  --input results/mass_spectrum_fine_scan.json --output results/mass_spectrum_reference.json
```

The bin config has one `mass_edges = ...` entry in GeV. `--scan` doubles FF
panels, doubles the grid twice, doubles span at fixed finest spacing, then
refines hard slivers/tolerances. It adds both fragmentation-scale and
starting-scale variations with common bin correlations. Each point records its
actual numerical controls. `--dump` exports the retained FF and hard source
arrays for an independent coarse audit. Saved checks hash both configurations,
the source, executable and attached scan.

Tests cover Born ownership on an exact M edge, empty above-M bins, every
additive component's closure, a direct real-tail fixed-order calculation,
unchanged shared FF quadrature, signed cell integration and invalid partitions.
Python independently constructs causal FFs and hard source cells, assigns
literal mass-cell overlaps, and checks the continuum O(alpha) bins with the
leading-power massive spectrum and allowed-azimuth geometry. It does not
independently reconstruct every resummed finite-geometry bin; that remains
covered by the inherited fixed-fraction and total-rate checks plus the new
per-bin refinements. See the [benchmark](../benchmarks/mass_spectrum_reference.md).
