# leptonff

**Perturbative timelike QED fragmentation and matched bare-lepton W/Z radiation.**
A C++17 numerical library with C, Fortran and Python interfaces. Its core
provides distribution-aware boundaries, kernels and evolution; matched modules
add spin-correlated finite decay responses. QCD generators connect through
optional adapters.

## Install

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build --parallel 4
cmake --install build --prefix "$PWD/install"
```

Use `-DBUILD_SHARED_LIBS=ON` for Python numerical bindings,
`-DLEPTONFF_BUILD_FORTRAN=ON` for a compiled Fortran module, or
`-DLEPTONFF_BUILD_MATCHED=OFF` for a standalone fragmentation/evolution package.
See the [installation guide](docs/manual/installation.md) for components and
client linking. Research reference programs are opt-in.

## First calculation

```cpp
#include <leptonff/core.hpp>
using namespace leptonff;
const auto ff = fragmentation_order1(
    Species::muon_minus, Species::muon_minus, QEDParameters{});
const auto momentum = ff.moment(2);
const auto endpoint_bin = ff.integral(.9L, 1);
const auto prediction = ff.convolve([](Real x) { return x*(1-x); }, .5L);
```

Link `leptonff::core` after `find_package(leptonff 0.26 CONFIG REQUIRED COMPONENTS core)`.
Matched W/Z clients include `<leptonff/matched.hpp>` and link `leptonff::leptonff`.
Complete examples cover [fragmentation](examples/fragmentation),
[convolution](examples/convolution), [matched W/Z](examples/cpp),
[C](examples/c_api), [Fortran](examples/fortran) and [Python](examples/python).

## Scientific documentation and checks

Start with the [scientific manual](docs/manual/index.md), especially the
[physics conventions and accuracy contract](docs/manual/physics.md),
[tutorials](docs/manual/tutorials.md), [interfaces](docs/manual/interfaces.md)
and [validation guide](docs/manual/validation.md).
Small [canonical benchmarks](benchmarks/canonical) ship with the source.

```sh
cmake -S . -B build -DBUILD_TESTING=ON -DLEPTONFF_BUILD_REFERENCE_PROGRAMS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --parallel 4
```

The optional MCFM 10.3.1 adapter supplies NLO QCD production and PHOTOS
comparisons. [Signed production archives and replay](docs/packaging.md) support
repeated QED/cut/bin studies; current MCFM exports contain central QCD weights.
Native DYTurbo integration and portable interpolated response tables remain
planned. No DYTurbo NNLO+NNLL result is available yet.

The photonic NLL decay prescription is explicitly unitary-conditioned. The Z
hard response supports parity-even, charge-symmetric measurements. Complete
off-shell EW and mixed QCD/EW corrections are outside scope. Numerical
validation does not by itself certify phenomenological bin precision.

## Contributing, citation and licence

See [CONTRIBUTING.md](CONTRIBUTING.md), [CITATION.cff](CITATION.cff) and
[CHANGELOG.md](CHANGELOG.md). Original code uses the [MIT licence](LICENSE);
the optional MCFM adapter retains [GPL-3.0-or-later terms](adapters/mcfm/LICENSE).
[Third-party notices](THIRD_PARTY_NOTICES.md) describe the licence boundaries.
Earlier studies remain in [development history](DEVELOPMENT_HISTORY.md).
