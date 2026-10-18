# Architecture and extension points

The package has three scientific layers:

1. **Core:** distributions, boundaries, kernels, coupling/flavour matching and
   evolution. It has no dependency on W/Z geometry, generator globals or PDFs.
2. **Matched processes:** finite W/Z hard terms, spin responses, measurements
   and explicitly conditioned photonic decay evolution.
3. **Applications:** signed production archives/replay, generator adapters,
   comparison scripts and research reference programs.

The core has a standalone exported library. The full library incorporates the
same core object files plus the matched/application functionality, retaining a
self-contained archive for existing host relinkers. Build options select the
needed functionality without introducing duplicate core compilation.

The source implementation of the special functions is independent of decay.
The x-space FF-to-decay validation helper resides in the matched layer. Its
header explicitly includes both evolution and decay geometry.

## Prepared physics and compiled measurements

`DecayCalculator` owns immutable physics. `MeasurementPlan` validates cuts and
bin definitions once and can be reused with compatible calculators. Its eta
selections share fraction actions. Converted-photon measurements share only
when the eta acceptance of every relevant parent direction agrees.

`FractionPairWorkspace` prepares exact atomic/cell cut stencils once per set of
bounds and grid. Several signed table channels reuse their searches, logarithms
and partial-cell factors. Scratch state is local to each evaluation; immutable
physics and measurements remain safe to share between workers.

The projected hard coefficient uses an analytic primitive shared with the
hard-table projection. Endpoint distributions, perturbative truncation and
the physical matching formula are unchanged.

## Exact-node response tables

`DecayResponseTable` builds responses concurrently at explicit parent momenta
and radiation seeds. It saves and reloads all five tensor coefficients per
rate/bin. Loading checks a complete configuration identity: software version,
numeric precision, evolution inputs, thresholds, grids, cutoffs, sampling,
cuts, bins, angular nodes and Cartesian frame convention. A footer and checksum
detect interrupted writes and accidental corruption. The checksum is not an
authenticity signature. Existing output files are not overwritten.

The table introduces no interpolation. Its indexed nodes may be contracted
with new compatible production tensors and PDF/scale weights. Callers must use
the stored parent momentum and frame, and separately establish angular and
hard-sampling convergence. Quadrature nodes must not be passed to the Monte
Carlo group estimator to manufacture statistical errors.

Response interpolation with measured errors and general evolution-table
serialization remain future work. Exact-node persistence provides the reusable
foundation without assigning an unverified interpolation accuracy.

## Replay and statistics

Replay batches complete signed groups and evaluates them on independent workers.
Results enter the estimator in archive order, preserving the arithmetic order
across worker counts and batch sizes. Every host block retains its own trial
denominator. Production weights reuse each computed model response.

The default full covariance is unchanged. Per-weight blocks retain all requested
bins/components/models within each production weight; diagonal mode retains
only variances. The C++ estimator additionally supports selected channel pairs.
Absent covariance entries are uncomputed, not zero. Details and commands are in
[the performance guide](performance.md).

## Generator integration

MCFM's adapter continues to use the full self-contained archive. The generic
matched interface accepts compatible Born points or production tensors, enabling
other hosts without a dependency in the core. DYTurbo still requires a native
adapter and a validated contraction of its production information with the
decay response. The package refactor does not supply an NNLO+NNLL result.

The design takes useful patterns from
[EKO](https://eko.readthedocs.io/en/latest/code/IO.html),
[HELL-x](https://www.roma1.infn.it/~bonvini/hell/),
[HOPPET](https://github.com/hoppet-code/hoppet) and
[Yadism](https://yadism.readthedocs.io/en/latest/theory/coeff-funcs.html).
These are design references, not added runtime dependencies.
