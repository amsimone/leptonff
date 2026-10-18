# leptonff

Lepton fragmentation and matched final-state QED radiation for bare-lepton W/Z
observables. A C++17 library with C/Fortran interfaces, spin-correlated decay
responses, and reusable signed QCD production records.

The research implementation includes photonic LL/NLL and coupled leptonic NLL
with flavour thresholds. An optional MCFM 10.3.1 adapter supplies NLO QCD W±/Z
production and PHOTOS comparisons. Native DYTurbo integration is planned;
**no DYTurbo NNLO+NNLL result is available yet**.

## Build, test and install

Requires CMake 3.16+, a C++17 compiler and, for tests, a C compiler and Python
3.10+. A Fortran compiler enables the installed Fortran example test. Core tests
require no generator, PDF installation or third-party Python packages.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --parallel 4
cmake --install build --prefix "$PWD/install"
```

Use `-DBUILD_SHARED_LIBS=ON` for a shared library. CMake consumers link
`leptonff::leptonff` after `find_package(leptonff 0.25 CONFIG REQUIRED)`.
See the [C++](examples/cpp), [C](examples/c_api) and [Fortran](examples/fortran)
examples. Native C++ clients need a compatible compiler ABI; the C interface
uses doubles at the ABI boundary.

## Reuse QCD production

Export pre-FSR signed production contributions with the MCFM adapter, then
change QED models, cuts and bins without rerunning QCD. Replay preserves
integration groups and cross-bin/cross-model covariance:

```sh
install/bin/leptonff-replay /path/to/production/manifest.json \
  --output predictions.json --subdivide 4 \
  --model photonic_ll --model photonic_nll --model leptonic_nll
```

See the [package and replay guide](docs/packaging.md) for export commands and
supported variations. Current MCFM archives contain central QCD weights only.
More bins require sufficient statistics and separate numerical checks.

## Documentation

- [Production/decay API and archive contract](docs/production_interface.md)
- [Photonic NLL conventions and matching](docs/photonic_nll.md)
- [Physics scope and accuracy audit](docs/qed_accuracy_audit.md)
- [Package and replay validation](benchmarks/package_reference.md)
- [DYTurbo status](docs/dyturbo_integration.md)
- [Development history](DEVELOPMENT_HISTORY.md) and [plan](plan.md)

This is research software. Implementation validation does not certify per-bin
phenomenological precision. The Z hard response supports parity-even,
charge-symmetric measurements. The decay prescription does not provide complete
off-shell EW or mixed QCD–EW corrections. Historical reports reference local
`results/` artifacts; large production datasets and external generators are
excluded from Git.

## Contributing, citation and licence

See [CONTRIBUTING.md](CONTRIBUTING.md), [CITATION.cff](CITATION.cff) and
[SECURITY.md](SECURITY.md). Original library code uses the [MIT licence](LICENSE).
The optional MCFM adapter retains [GPL-3.0-or-later terms](adapters/mcfm/LICENSE).
See [third-party notices](THIRD_PARTY_NOTICES.md) for details.
