# DYTurbo 1.4.2 interface

The native **QCD production client** is implemented. It calls DYTurbo's library
integration terms, writes machine-readable cross sections in fb, and exposes
NLO/NNLO with quadrature or VEGAS and NLO+NLL/NNLO+NNLL resummation. The original
executable is retained as an independent steering/reference check.

**The leptonff radiative decay is not applied by this production client.** These are
QCD production benchmarks with native Born-lepton decays, not bare-muon QED
predictions. The generic leptonff spin response and exact-node tables are
available. A separate [differential W decay client](dyturbo_w_decay.md) now
applies them to native NLO+NLL RES/CT/V+jet points, with PHOTOS and QCD export.
Its pilot validation does not certify a high-precision W shape result.

A separate [inclusive dimuon-mass example](dyturbo_qed_mass.md) now convolves
NLO+NLL QCD production with the order-alpha and matched photonic NLL mass
responses. It uses a refined midpoint mass-migration approximation and no
laboratory lepton cuts; it does not require the general angular interface.

## Build and run

No generator, PDF grid, generated executable or dependency is vendored. The
optional bootstrap uses checksum-pinned official release archives and puts
sources, installed libraries, logs and a build manifest in an ignored directory:

```sh
python3 adapters/dyturbo/build.py --jobs 4
python3 adapters/dyturbo/run.py --output results/dyturbo-z \
  --pdf-data /path/to/LHAPDF/grids --calls 100000 --accuracy .001
```

Requirements: Python >=3.10, curl, make, CMake, C/C++ and GNU Fortran compilers.
On macOS the bootstrap uses Apple Clang with GNU Fortran, disables OpenMP and
the platform-incompatible TLS optimisation, locates the Fortran runtime, and
uses ordinary optimisation without fast-math. The serial build is intentional;
parallel generator execution has not been certified. Run separate processes
in separate output directories for independent jobs.

The bootstrap builds LHAPDF 6.5.4 with the legacy C++ interface still required
by this release, Cuba 4.2.2, cubature 1.0.3 and HELLx 3.0. HELLx is a source/link
dependency; its presence does not mean small-x resummation is enabled. Two
HELLx template fixes follow the upstream DYTurbo installer; two DYTurbo VLA
initialisers are replaced by equivalent explicit zeroing for Clang. All edits
are isolated and their resulting hashes are recorded. Dependency archive
extraction rejects traversal and escaping links and preserves autotools
release timestamps. Cuba is built serially to avoid archive-member races.

For an existing DYTurbo installation, build only the client:

```sh
cmake -S adapters/dyturbo -B build-dyturbo-interface \
  -DDYTURBO_ROOT=/path/to/dyturbo/prefix
cmake --build build-dyturbo-interface --parallel 4
```

The CMake project checks the backend version and uses that prefix explicitly.
This adapter stays outside the standalone library's dependency graph.

## Supported production matrix

| Mode | Fixed-order matching | QCD logarithms | Integration |
| --- | --- | --- | --- |
| `nlo-quad` | NLO | none | quadrature |
| `nlo-vegas` | NLO | none | VEGAS |
| `nnlo-quad` | NNLO | none | quadrature, inclusive leptons |
| `nnlo-vegas` | NNLO | none | VEGAS |
| `nlo-nll` | NLO | NLL | quadrature |
| `nnlo-nnll` | NNLO | NNLL | quadrature, inclusive leptons |

With no `--mode`, the runner executes all six. Repeat `--mode` to select a
subset. Processes are `z`, `wp`, `wm`; use e.g. `--process wp --mass 60 100`.
Mass, boson transverse momentum and rapidity bounds can be set explicitly.
All runs use 13 TeV pp collisions. The PDF set/member are explicit. The default
NNPDF31 NNLO member is shared across orders to compare integrators; these are
not separate order-consistent PDF-fit comparisons.

The runner fixes the Sudakov, hard, evolution, exponentiation, counterterm and
V+jet orders explicitly and disables `orderlha`, so the selected PDF cannot
silently change the resummation order. DYTurbo 1.4.2 calls the included hard
matching convention “primed”; `primed=true` retains the hard and subtraction
coefficients at the requested fixed order. The labels here describe the
logarithmic towers and fixed-order matching separately. The resolved orders
are exported and checked. DYTurbo's own QED/QCD-QED resummation is disabled;
it is separate from final-state lepton fragmentation.

`--fiducial` enables native Born-lepton cuts from the release default card.
DYTurbo 1.4.2 explicitly rejects NNLO quadrature with lepton cuts because its
analytic NNLO V+jet contribution is angular-inclusive. The runner rejects that
request before starting any run. `--hybrid` explicitly permits quadrature for
the available pieces and VEGAS for NNLO V+jet. Per-term method labels preserve
that distinction. A hybrid result must not be called fully quadrature.

## Output and verification

The [local native-backend benchmark](../benchmarks/dyturbo_native_1_4_2.json)
records all six modes for Z and W+, refinement checks and the isolated build.
For Z at 13 TeV with 66 < m < 116 GeV and the shared NNLO PDF, refined NLO
quadrature gives 1884.404 +/- 0.169 pb in 2.88 s; VEGAS gives 1882.587 +/-
6.28 pb in 6.31 s. Refined NNLO quadrature gives 1909.586 +/- 0.747 pb in
221 s; the 300,000-call VEGAS check gives 2159.299 +/- 155 pb in 268 s.
These errors are backend estimates and the timings include different concurrent
loads. Analytic production does not eliminate numerical integration, and the
NNLO VEGAS result is still a low-precision check. All six Z client results,
both refined NLO modes, and two runs of the independently bootstrapped backend
reproduce the original executable to the checker's numerical tolerance.

Each mode has its input card, backend log, raw `native.json` and checked
`result.json`; `matrix.json` collects the runs. Cross sections are integrated
values in fb, with signed terms retained. Resolved settings include the QCD
orders, integration dimensions, PDF/member, seed, EW masses/widths, scales,
subtraction cutoff and resummation prescription. The build manifest records
source archives, changes and compiler commands. New client builds also record
backend-library/default-card hashes; explicit PDF paths record the PDF info
and member-file hashes. Output directories and result files are not overwritten.

The runner checks term sums, perturbative orders and integration methods. It
uses DYTurbo's own combination of term uncertainties; quadrature nodes are
never fed into an event Monte Carlo covariance estimator. Finite-budget error
estimates are not a certified accuracy guarantee. In particular, NNLO VEGAS
has large real/virtual cancellations and needs substantially more statistics
than an execution smoke test. Quadrature errors require refinement; VEGAS
requires independent seeds; both require subtraction-cutoff scans.

To compare the client against the original executable with identical physics
settings and random seeds:

```sh
python3 adapters/dyturbo/check_reference.py \
  --results results/dyturbo-z \
  --executable external/dyturbo-native/source/dyturbo-1.4.2/local/bin/dyturbo
```

The checker reads the original full-precision text result, checks every bin's
value and uncertainty, and writes closure reports. This establishes native
client closure in QCD; it does not validate a QED coupling. Lightweight CI
checks steering, signed sums, method/order guards and safe source preparation
without downloading or simulating DYTurbo. Real-backend results must be
identified separately from those contract tests.

## Remaining radiative coupling and precision work

The [DYTurbo paper](https://arxiv.org/abs/1910.07049) describes the production /
angular-acceptance factorisation underlying the fast integration. Source 1.4.2
has distinct angular moments for Born/resummed/counterterm terms, differential
LO V+jet integration, and angular-integrated analytic NNLO V+jet. The last
path cannot supply arbitrary W lepton cuts by multiplying a scalar correction
onto a finished histogram.

The W differential client maps the native charge ordering to leptonff's
Cartesian spin response, replaces the measurement in RES/CT/V+jet, and passes
zero-QED and production-replay closure. The next fast route should connect
native angular moments directly to the tensor response, with correlated
finite-term integration. General Z fiducial and NNLO coupling remain open.
Mapped subtraction points need common radiation streams. No term may divide
by a Born weight or clip signed contributions. For NNLO fiducial spectra the
initial supported route will need the differential VEGAS/hybrid backend;
full deterministic NNLO fiducial integration needs additional angular-resolved
production information. NNLO quadrature remains useful for inclusive checks
and rotationally invariant measurements with a separately validated decay
factorisation.

Then test photonic order-alpha and LL/NLL limits and produce bare-lepton W pT /
eta and Z mass spectra against a shared PHOTOS reference. Production caching,
PDF/scale-weight vectors, interpolation error controls and a high-statistics
phenomenology campaign remain separate work. The W response is factorised /
pole-type and the neutral response is parity-even for charge-symmetric
measurements. This does not provide full off-shell EW or complete mixed QCD-EW
corrections.
