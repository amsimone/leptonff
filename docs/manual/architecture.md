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

## Reuse and future tables

Current evolution and decay objects can be constructed once and reused in
memory. A portable future table must identify its schema, software version,
physical settings, basis, normalization, valid kinematic domain, numerical
controls and interpolation error. Evolution tables and fiducial decay-response
tables are different objects: the latter also depend on process, cuts, bins,
spin basis and supported variations.

Portable table serialization and response interpolation are not implemented
in this version. Their acceptance criterion is agreement with the direct
calculation under independent grid/kinematic refinement, together with a measured
benefit in repeated evaluation time.

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
