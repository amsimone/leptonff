# Contributing

Open an issue for bugs, physics questions or proposed changes. For numerical
issues include the version/commit, compiler, build options, configuration,
random seeds, expected result and the smallest reproducer available. Separate
statistical errors from quadrature, grid and theory uncertainties.

Build and run the standalone checks before submitting a pull request:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DLEPTONFF_BUILD_REFERENCE_PROGRAMS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --parallel 4
```

Tests include independent numerical references, zero-coupling/matching limits,
signed-event covariance, replay and relocated package consumers. Optional host
campaigns require separate installations and are not run by public CI. Plotting
additionally needs NumPy and Matplotlib; core tests use Python's standard library.

Start with the [scientific manual](docs/manual/index.md) and its interface policy.
For physics changes document approximations, coupling/normalisation conventions,
source references and relevant limits. Add meaningful independent or regression
checks where appropriate. Preserve signed subtraction groups and common radiation
seeds when changing production/replay interfaces.

Keep production records, build outputs, external generator sources, credentials
and machine-specific configuration out of commits. Historical `results/` links
refer to separately retained local artifacts. Summarise changes and validation
in the pull request; retain third-party notices. Contributions use the licence
of the files changed: MIT by default, GPL-3.0-or-later in `adapters/mcfm/`.
