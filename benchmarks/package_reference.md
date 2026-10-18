# Version 0.25 package and replay validation

The generator-independent library, C/Fortran ABI, configurable laboratory
bins, signed production archive and multi-model replay are implemented.
Usage and limitations are in [the package guide](../docs/packaging.md).
This benchmark validates implementation equivalence, not phenomenological
precision or a native DYTurbo coupling.

## Package checks

The full regression run passed 69/69 tests before the last diagnostic and
installation changes. After those changes, all nine affected checks passed,
including two newly registered tests: the five-model replay CLI and safe
source preparation. The configured suite now contains 71 tests; a single
full 71-test run is not claimed. Evidence:

- [Full regression log](../results/package_ctest_v025.log).
- [Final affected checks](../results/package_final_checks_v025.log).
- [Shared-library installation check](../results/package_shared_install_v025.log).

Both static and shared installations were moved to a different directory.
Independent C, C++ and Fortran consumers built and evaluated a decay using
the installed package; the installed replay executable also started there.
The archive/covariance tests cover signed subtraction groups, omitted zero
trials, incomplete files and resumed iterations with different call counts.
The source preparer was tested on synthetic archives only.

Subsequent repository-preparation checks passed the complete 71-test suite
locally with GCC (static) and Apple Clang (shared), including the installed
consumers. See [the GCC log](../results/github_gcc_tests.log) and
[the Clang log](../results/github_clang_tests.log). A C++17-only test syntax
cleanup was checked separately afterward. Linux GCC static/shared and macOS
Clang jobs are configured in GitHub Actions; a hosted CI result is not yet
available. Historical log and dataset links refer to retained local artifacts,
not files bundled in the source repository.

## MCFM online/offline closure

Each comparison checks all five analytic channels across the total and 17
observable bins (90 channels), including the full covariance matrix. These
are reduced-call smoke runs with identical records and radiation seeds.

| Process | Maximum mean difference [pb] | Maximum covariance difference [pb²] |
| --- | ---: | ---: |
| W+ | 1.82e-12 | 2.18e-11 |
| W−, fixed calls | 2.61e-12 | 1.75e-10 |
| Z | 2.72e-13 | 9.09e-13 |
| W−, resumed/variable calls | 2.05e-12 | 1.46e-10 |

Reports: [W+](../results/package_replay_wplus_checked_v025.json),
[W−](../results/package_replay_wminus_checked_v025.json),
[Z](../results/package_replay_z_checked_v025.json),
[resumed W−](../results/package_replay_online_resumed_checked_v025.json).
The resumed run exercises 24 contribution/iteration blocks and varying
trial denominators. An earlier W− pilot stopped at the old stage-one limit;
its failed manifest is retained and excluded from this closure evidence.

Fourfold finite-bin subdivision was exercised offline (62 W bins; the
configured Z definition has 65). Infinite overflow bins remain unsplit.
The all-model CLI test uses synthetic records; it is not a five-model
high-statistics QCD campaign. Existing v0.23 comparison figures are unchanged.

## Remaining work

- Obtain and inspect the requested DYTurbo 1.4.2 source, implement its native
  acceptance/tensor coupling and demonstrate QED-off and matched closure.
  See [the source and integration status](../docs/dyturbo_integration.md).
- Establish a differential precision budget and launch suitable independent
  production samples. Fourfold bin subdivision does not establish tenfold
  smaller errors. Host convergence warnings remain in the run manifests.
- Validate decay quadrature, hard-sampling and any future response-table
  interpolation separately from QCD integration covariance.
- Export explicit scale/PDF weight vectors for offline variations. The
  current MCFM exporter stores central weights only.

There is no native DYTurbo NNLO+NNLL result, certified response grid or new
high-statistics phenomenology plot in this increment.
