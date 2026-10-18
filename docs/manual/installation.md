# Installation and components

The native library requires CMake 3.16+ and a C++17 compiler. Python 3.10+
enables reference checks and the optional Python bindings. No generator, PDF
installation or third-party numerical library is required for the core.

## Native installation

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build --parallel 4
cmake --install build --prefix "$PWD/install"
```

This builds the core, matched W/Z library, introductory examples and replay.
Research reference programs are opt-in. Useful options are:

| CMake option | Default | Effect |
|---|---|---|
| `BUILD_SHARED_LIBS` | `OFF` | Shared libraries; required for Python numerical bindings |
| `LEPTONFF_BUILD_MATCHED` | `ON` | W/Z matching, geometry and production records |
| `LEPTONFF_BUILD_TOOLS` | `ON` | Installed replay executable and Python replay command |
| `LEPTONFF_BUILD_EXAMPLES` | `ON` | Introductory fragmentation, convolution and matched examples |
| `LEPTONFF_BUILD_REFERENCE_PROGRAMS` | `OFF` | Research executables and their reference tests |
| `LEPTONFF_BUILD_FORTRAN` | `OFF` | Compiled Fortran module and exported target |
| `BUILD_TESTING` | `ON` | Unit tests and available interface/installation checks |
| `LEPTONFF_BUILD_DOCS` | `OFF` | HTML manual; requires Sphinx and MyST |

For a small fragmentation/evolution installation, add
`-DLEPTONFF_BUILD_MATCHED=OFF -DLEPTONFF_BUILD_EXAMPLES=OFF`.
This also disables replay and the decay reference programs.

## C++ clients

```cmake
find_package(leptonff 0.27 CONFIG REQUIRED COMPONENTS core)
target_link_libraries(my_program PRIVATE leptonff::core)
```

Include `<leptonff/core.hpp>`. Matched clients request `COMPONENTS matched`,
include `<leptonff/matched.hpp>` and link `leptonff::leptonff`.
The complete library includes the core, preserving the self-contained
`libleptonff.a` archive used by existing generator relinkers. Choose one of the
two library targets for a given executable; linking both is unnecessary.

`pkg-config --cflags --libs leptonff` supports ordinary C++ Makefiles. Link with
a C++ linker, including when the calling source uses the C ABI.

Native C++ uses `long double`; the installed target checks its mantissa ABI.
C, Fortran and Python use doubles at the ABI boundary. The installed documentation
includes the canonical analytic benchmark data.

## Fortran

Configure with `-DLEPTONFF_BUILD_FORTRAN=ON` and, when necessary,
`-DCMAKE_Fortran_COMPILER=gfortran`. Clients use `use leptonff` and link
`leptonff::fortran`. Compiled module files are compiler-specific. The installed
source is also available through the CMake variable `leptonff_FORTRAN_SOURCE`;
the standalone Fortran example demonstrates compiling it with the client's own
compiler when a compiled module was not installed.

## Python

Build and install a shared library, then install the numerical bindings:

```sh
python3 -m pip install .
export LEPTONFF_LIBRARY="$PWD/install/lib/libleptonff.so"
python3 examples/python/fragmentation.py
```

On macOS use `libleptonff.dylib`. A core-only application can select
`libleptonff_core.so` or `.dylib`. The Python package uses the installed native
library through its versioned C ABI; it requires no NumPy dependency and does
not bundle or download a generator.

## Manual and API reference

```sh
python3 -m pip install -r docs/manual/requirements.txt
sphinx-build -W --keep-going -b html docs/manual build/manual
```

Open `build/manual/index.html`. Alternatively enable `LEPTONFF_BUILD_DOCS` in
CMake; CMake installations then include the generated HTML manual. With Doxygen installed, `doxygen Doxyfile` generates the C/C++ symbol
reference under `build/api/html`.
