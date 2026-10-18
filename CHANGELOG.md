# Changelog

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

Native DYTurbo integration and a certified high-statistics phenomenology campaign
remain open. See [the history](DEVELOPMENT_HISTORY.md) for earlier milestones.
