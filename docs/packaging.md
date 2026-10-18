# Installing and using leptonff

Version 0.25 separates the generator-independent decay calculation from the
MCFM adapter. The package contains a C++17 library, a versioned C interface,
a standalone production replay program and a checksum-verifying Python
front end. MCFM, DYTurbo and PHOTOS are not library dependencies.
See [the validation benchmark](../benchmarks/package_reference.md) for the
static/shared installation checks and online/offline closure evidence.

## Build and install

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON
cmake --build build-release --parallel 4
ctest --test-dir build-release --output-on-failure --parallel 4
cmake --install build-release --prefix "$PWD/install"
```

Static builds remain available with `BUILD_SHARED_LIBS=OFF` (the default).
The existing MCFM relinker needs a static core build using the host's
compatible GCC/libstdc++ ABI. It is deliberately separate from installation.

Consumers use `find_package(leptonff 0.25 CONFIG REQUIRED)` and link
`leptonff::leptonff`. See `examples/cpp`, `examples/c_api` and
`examples/fortran`. The installation is relocatable; the examples are tested
against a moved installation directory, without source/build include paths.

The native C++ interface uses `long double` and the compiler's C++ ABI.
An exported mantissa check rejects a different long-double format. Use a
compatible compiler/standard library, or the double-based C interface in
`leptonff/c_api.h`. Its opaque contexts, explicit array capacities and
exception containment are suitable for C and `ISO_C_BINDING` callers.
`leptonff_decay_config_default` initialises the ABI version and structure size.
Contexts are immutable; evaluation and errors are local to each calling
thread. C diagnostics expose the photonic source-normalisation convention.

## Generate QCD once, replay decay studies

For the existing MCFM installation, build the optional adapter as before,
with `--iteration-audit`. Then export production without running any QED:

```sh
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-portable-v025/mcfm_leptonff \
  --output results/my_wplus_production --threads 1 \
  --export-production --production-only
```

The call budget, PDF/scales and technical subtraction settings still need to
be chosen and validated for the desired precision. This command is not an
automatically high-statistics or precision-certified setup. Repeat for W-
and Z; these are separate production datasets.

With the installed bin directory on PATH, replay the completed manifest:

```sh
leptonff-replay results/my_wplus_production/manifest.json \
  --output results/my_wplus_qed.json --subdivide 4 \
  --model photonic_ll --model photonic_nll --model leptonic_nll
```

In the source checkout the same front end is
`python3 python/replay_production.py ... --engine build-host/leptonff_replay`.
The Python front end verifies the archive SHA256, original configuration,
pre-FSR/full-angle contract and generated Z mass coverage. It preserves the
host's integration warnings and records the replay executable and arguments.
It also rejects production inputs changed during evaluation.

`leptonff_replay` is the lower-level C++ program and accepts a `.lff` file
directly. That route validates the record structure and completion marker;
the producer remains responsible for its physics/coverage declarations.

All models in one invocation consume the SAME signed records and radiation
seeds. The output contains cross-bin AND cross-model covariance, not just
separate errors. No rerun of QCD or PHOTOS is involved. PHOTOS remains an
optional external benchmark through the MCFM adapter.

## Models and observables

Available replay models are `born`, `photonic_ll`, `photonic_nll`,
`leptonic_nll` and `qed_hvp`. The last uses the declared pointlike/HVP input
prescription with strong FF evolution off. Numerical/physical content is
described in the accuracy audit and photonic NLL contract; extra fermion
content is not an additive cross section.

The photonic NLL source is explicitly conditioned to the unitary decay
normalisation. Its removed inclusive factor, source first moment and domain
bounds are exported per model. It is not an absolute NNLO partial width.
The `order2` channel is FF-only in photonic modes (NLL adds boundary/kernel
terms), and the full retained matching coefficient in coupled modes.

The output is in pb and pb^2. It includes total-rate and bin means, full
integration covariance, actual model/grid/cut configurations and every
iteration's call count. Equal weights are used across production iterations
within each contribution. Resumed stages and differing call counts retain
their own zero-trial denominators. Adaptive stopping and rare signed tails
still require independent-seed diagnostics; covariance is not a theorem of
precision or a theory-uncertainty band.

Fourfold subdivision applies to finite bins, leaving infinite overflow
unsplit. W has 62 bins instead of 17; Z has 65. Thus the finite W pT and eta,
and Z mass distributions each gain exactly four times as many bins.
This does NOT reduce relative statistical errors: those require sufficient
production and decay sampling, and numerical refinement separately.

An explicit binning file has one observable per line:

```text
pt 25 26.25 27.5 28.75 30 35 40 45 60 inf
abseta 0 0.15 0.30 0.60 1.2 1.8 2.4
```

Pass it with `--binning edges.conf`; do not also pass `--subdivide`.
Supported names are W `pt`, `mt`, `abseta`, and Z `mass`, `minpt`, `maxeta`.
Only requested observables are calculated. Edges must partition their full
fiducial support, including pT/mT overflow; nested strict cuts assign internal
boundaries once. Changing fiducial cuts requires compatible explicit edges.
Options include `--pt-min`, `--eta-max`, `--mass-min`, `--mass-max`,
`--nu-pt-min`, `--mt-min`, `--mu-frag`, `--grid-bins`, `--hard-samples` and
the two finite-geometry cutoffs. Wider selections require generated support.

The current MCFM exporter stores CENTRAL QCD weights only. Its records do
not reconstruct arbitrary PDF, QCD-scale, collider-energy or mass changes.
The format supports multiple named weight vectors for future exporters.
Existing total/spectrum CSV files are not reusable production archives.

For other generators and direct spin-tensor calculations, see
`production_interface.md`. DYTurbo-specific status is in
`dyturbo_integration.md`.
