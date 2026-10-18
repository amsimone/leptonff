# Supported interfaces

## C++

| Entry point | Objects and operations | Library target |
|---|---|---|
| `leptonff/core.hpp` | Distributions, fixed-order fragmentation, LL/coupled evolution, physical inputs | `leptonff::core` |
| `leptonff/matched.hpp` | Decay presets/calculators, compiled measurements, response tables, Drell–Yan rotations | `leptonff::leptonff` |
| `leptonff/production_records.hpp` | Signed groups, archive I/O, replay covariance | `leptonff::leptonff` |

The component headers included by these entry points remain available.
Named physics objects and documented operations form the supported interface.
Raw lattice arrays, validation helpers and historical spectrum constructors
are advanced research interfaces and may change between minor releases.
The C ABI has its own explicit version and layout checks. Record the runtime
`leptonff::version()` (C: `leptonff_version()`, Python: `context.native_version`)
alongside numerical and theory settings for reproducibility.

Rebuild C++ clients when moving to 0.27: matched implementation layouts changed.
The full shared library now uses a minor-version SONAME, matching the CMake
package's minor-version compatibility boundary. Existing C structure layouts
remain ABI-v1; newly linked clients use the new full-library SONAME. The core's
shared-library ABI is unchanged.

`EvolutionTheory` contains the alpha scheme/reference, masses, scales, running
mode and optional flavours. `MomentEvolutionOptions` and `QEDXSpaceGrid`
independently control numerical refinement. Existing explicit backend parameter
structs remain available to advanced users.

`Distribution::integral`, `moment`, `apply` and `convolve` give complete
distribution actions. Neither `resolved_density` nor a table of interior values
alone represents the endpoint distribution.

## C

`leptonff/fragmentation_c.h` is usable with the standalone core. Initialize a
`leptonff_fragmentation_config` using its default function, edit named fields,
create an opaque context, and call moment/integral/convolution operations.
Destroy the context when finished. Species IDs are 0 (negative muon),
1 (positive muon), 2 (photon).

The default model is full first-order one-family fragmentation. Coupled NLL
accepts continuum moments with `N>1`. Set `grid_bins>0` at creation to prepare
coupled x-space bin actions; those require a positive lower bound. Fixed-order
convolutions accept a callback with user data.

`leptonff/c_api.h` additionally exposes the matched decay interface. New clients
use `LEPTONFF_PROCESS_W_PLUS`, `LEPTONFF_PROCESS_W_MINUS` and
`LEPTONFF_PROCESS_Z`, whose values are PDG boson IDs. Legacy ABI-v1 process
values 1, 6 and 31 remain accepted.

For compatibility, the ABI-v1 decay preset has one alpha field: fixed alpha in
photonic modes and on-shell alpha in leptonic/HVP presets. The fragmentation
configuration instead makes the scheme and reference scale explicit.
`leptonff_decay_create_with_theory` adds explicit scheme, reference/initial scale
and running order without changing the ABI-v1 preset layout. MSbar alpha must
be provided at the initial scale. Photonic evolution requires fixed alpha and
initial scale equal to the muon mass; the HVP anchor requires on-shell alpha.
When explicit theory is supplied, its alpha overrides the preset alpha field.

`leptonff_decay_with_observables` creates a view sharing the physics context;
only cuts and binning may change. Either handle may outlive the other.
`leptonff_decay_evaluate_batch` accepts contiguous `(count,4)` momentum arrays
and one seed per point. Its output is event-major, then rate/bin; each entry
contains the five existing prediction fields. Capacity counts predictions.
Evaluation errors are joined and returned before writing the caller's output.
String queries expose the resolved physics and ordered output labels.

Fallible operations return zero on success and nonzero on failure; retrieve
the thread-local message with `leptonff_last_error`. Exceptions never cross the
C boundary. Caller-provided callbacks must return finite values and must not
throw through C. Contexts are immutable; coordinate destruction with users of
the same handle.

Decay result fields have the following meaning:

| Field | Meaning |
|---|---|
| `born` | Born acceptance/response |
| `order1` | Complete implemented prediction through first order, including Born |
| `matched` | Combined fixed-order and resummed prescription |
| `order2` | Retained second-order term; a diagnostic, without a complete NNLO claim |
| `higher` | `matched - order1` |

These are related predictions and diagnostics. Apply the host production weight
when forming a cross section; do not sum the five output fields together.

## Fortran

The installed `leptonff` module declares C-interoperable fragmentation and decay
configuration/prediction types, named constants and lifecycle/evaluation calls.
The module centralizes declarations previously copied into the example.
Use `integer(c_size_t)` for counts and `real(c_double)` for ABI values.
The installed example checks `c_sizeof(config)` against the native size before
creating a context. No C++ `long double` layout is exposed through Fortran.
Batch momenta use `first(4,count)` / `second(4,count)`, contiguous by point;
seeds use `integer(c_int64_t)`. The module also declares explicit theory,
measurement-view and string-query entry points.

## Python

`Fragmentation` and `Decay` own native contexts and support `with` blocks and
explicit `close()`. Inputs are frozen configuration dataclasses; query methods
return Python numbers, prediction dataclasses or complex tensors. Errors from
the library become `ValueError`; using a closed handle raises `RuntimeError`.
Python callback exceptions are caught before returning through C, then raised
again in Python. The binding uses the same native numerical implementation.

`Decay.with_observables` shares prepared physics with a new measurement.
`evaluate_batch` returns ordinary spectra from one native call. Optional
`evaluate_array` accepts NumPy arrays and returns float64 shape
`(points, outputs, 5)`, with final axis `born,order1,matched,order2,higher`.
`output_labels` names the second axis. `resolved_configuration()` records the
actual inputs; `DecayTheory` supplies explicit coupling settings when needed.
NumPy is optional (`pip install 'leptonff[arrays]'`).

## Migration from 0.25

- Existing `leptonff::leptonff` and the self-contained full archive remain.
- `decay_xspace_ff_match` and `XSpaceFFDecayMatch` moved to
  `leptonff/ff_decay.hpp`; include it explicitly in validation clients.
- Reference executables retain names, and are built under `build/apps/` when
  enabled. Replay remains in the root build directory and installed `bin/`.
- Standalone examples request CMake package version 0.27.
- The historical low-level headers remain installed; the new manual defines
  the supported entry points and distinguishes advanced controls.
