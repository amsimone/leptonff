# Changelog

## Unreleased

- Add a pinned optional McMule strict O(alpha^2) photonic shape benchmark,
  inclusive NLO normalization, independently checked importance maps, six-Q
  logarithm fits and FKS/technical controls. Find finite moment remainders
  consistent with the retained NLL mass logs; keep full NNLO matching and
  resummed residual interpretation distinct. Document computation ingredients,
  external evidence, uncertainty limits and missing phenomenology inputs.

- Audit the complete normalized Z radiator through second order with independent
  continuum moments/cumulatives, mass/coupling scans and numerical controls.
  Add a pinned ROOT-free KKMC adapter with paired order weights and controls,
  a mass-only beyond-NLL matching variation, and separate uncertainty diagnostics.
  Report resolved external residuals without assigning a total theory error.
  Support independent multi-seed KKMC campaigns and covariance-preserving raw
  sample pooling, with seed consistency diagnostics and tenfold reference statistics.
  Add coupling/mass residual diagnostics with absolute fit covariance and
  truncation checks; expose the legacy CEEX support floor and guard sample pooling.

- Add leptonff-only fine-bin Z mass accuracy plots and a QED fragmentation-scale
  scan with consistent hard matching, first-order scale-cancellation tests,
  separate grid/quadrature controls and magnified NLL/LL differences.

- Add an optional ROOT-free ReneSANCe/SANC Born-QCD FSR reference, independent
  massive Z mass radiator with quadrature, W fixed-order normalization and
  cutoff checks, PHOTOS/PYTHIA decay comparisons, and original/additive-NLO
  matching plots. Extend the inclusive mass reduction to matched LL as well
  as NLL. Fine-bin shower/W statistical precision remains explicitly open.

- Demonstrate photonic order-alpha and matched NLL decay on common DYTurbo
  NLO+NLL QCD production for the angular-inclusive dimuon mass spectrum,
  with mass-migration and fragmentation-grid refinements and zero-alpha
  closure. Add a fine-bin conservative mass convolution and reusable PHOTOS
  transfer moments for a shared-production comparison.
- Implement differential DYTurbo NLO+NLL W decay coupling in RES/CT/V+jet,
  correlated per-bin model covariance, QCD production export and a runner.
  Zero-alpha/native/replay pilot checks pass; precision W shapes remain open.

- Add an isolated, checksum-pinned DYTurbo 1.4.2 build and native QCD production
  client, a six-mode NLO/NNLO/resummed integration runner, signed per-term JSON
  results and an original-executable closure checker. NNLO fiducial quadrature
  is explicitly guarded.
- Resolve the leptonic decay preset's observed-flavour threshold from the
  caller's muon mass, including C ABI promotion on extended-precision systems.

## 0.27.0

- Replace event-loop hard-kernel quadrature with a shared analytic primitive.
- Reuse exact signed fraction-pair cut stencils and compiled eta measurements.
- Add deterministic batch/parallel evaluation and bounded grouped replay.
- Add block, selected-pair and diagonal covariance, preserving full mode by default.
- Share prepared physics across C/Python measurement views; expose resolved
  settings, output labels and explicit matched alpha inputs.
- Add C/Fortran/Python batch interfaces and optional NumPy arrays.
- Persist exact-node spin responses with complete configuration identities,
  round-trip precision, truncation checks and a corruption checksum.
- Keep physical matching prescriptions and legacy C structure layouts intact.
  C++ clients should rebuild for the new minor version. The full shared library
  now uses a minor-version SONAME; the standalone core retains its existing ABI.

## 0.26.0

- Standalone fragmentation/evolution target `leptonff::core`; the full
  self-contained archive remains compatible with generator relinkers.
- Runtime native version reporting through C++, C/Fortran and Python.
- Explicit alpha scheme/reference inputs and separate theory/numerical settings.
- Endpoint-aware distribution bin integrals and coefficient convolutions.
- Versioned core C API, importable Python numerical bindings, and an installed
  Fortran module/source replacing duplicated declarations in examples.
- Scientific manual, three introductory C++ examples, canonical analytic
  benchmarks, pkg-config metadata and opt-in research executables.
- Splitting kernels no longer depend on the vector-decay source. The x-space
  decay validation helper moved to `leptonff/ff_decay.hpp`.
- Existing numerical physics and ABI-v1 decay semantics are retained. Portable
  response tables and native DYTurbo integration remain future work.

## 0.25.0

- Installable static/shared C++17 library, CMake target and versioned C interface,
  with C, C++ and Fortran examples.
- Generator-independent weighted-point and spin-tensor decay responses.
- Configurable bins and signed production archives with grouped covariance,
  including resumed iterations.
- Checksum-verified multi-model replay using shared production.
- MCFM production-only export and online/offline W±/Z closure.
- Repository licensing, contribution guidance and standalone CI.

A certified high-statistics fiducial phenomenology campaign remains open;
the native DYTurbo coupling is now demonstrated in the configurations above.
See [the history](DEVELOPMENT_HISTORY.md) for earlier milestones.
