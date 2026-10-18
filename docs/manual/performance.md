# Fast evaluation and repeated studies

Prepare a calculator and a `MeasurementPlan` once. Reuse them for every point;
the overload accepting a raw `LabSpectrum` remains available and compiles that
measurement on each call. Physics contexts contain no production PDFs or weights.

```cpp
#include <leptonff/matched.hpp>
using namespace leptonff;
auto settings = make_decay_configuration(SpinDecay::wplus,
    DecayPreset::photonic_nll, AlphaInput::on_shell(1/137.035999084L),
    .1056583755L, 80.385L);
DecayCalculator decay(settings);
MeasurementPlan measurement(false, LabSpectrum::configured(false,
    LabSpectrum::reference_binning(false, 4)));
auto result = decay.evaluate_batch(points, measurement, 4);
```

Here `points` is a vector of `BornDecayPoint`. Preserve a common radiation seed
for real/mapped counterevents. Batch evaluation returns functionals, not weighted
cross sections or integration errors.

## Offline replay

```sh
leptonff-replay manifest.json --output fine.json \
  --model photonic_nll --model leptonic_nll --subdivide 4 \
  --threads 4 --batch-groups 128 --covariance full
```

The recommended manifest front end continues to check production provenance.
The native `leptonff_replay` engine accepts the same evaluation options directly
for an archive. Thread count defaults to one. Batches are bounded to roughly
64 MiB of result vectors, in addition to the stored production records and
covariance. Worker exceptions stop the batch; no incomplete result is emitted.

| Covariance option | Retained information | JSON format |
|---|---|---|
| `full` (default) | All bins, models, components and production weights | Schema 1, dense `covariance_pb2` |
| `per-weight` | Full covariance within each production weight | Schema 2, `covariance_pairs` and `covariance_values_pb2` |
| `diagonal` | Channel variances only | Schema 2, diagonal pairs |

All formats include `variance_pb2`. Schema-2 pairs index `labels`; entries not
listed were not computed. Per-weight output cannot supply correlations between
PDF/scale variations, and diagonal output cannot supply bin/model difference
errors. Choose the layout for the intended analysis. Full storage is limited
to 8,192 channels; sparse layouts accept larger vectors, subject to an explicit
stored-entry memory guard. The C++ `CovarianceRequest` also accepts selected
cross-channel pairs; its output always includes every diagonal.

## Python batches and measurement views

```python
import math
from leptonff import Decay

with Decay() as physics:
    with physics.with_observables({"pt": [25, 30, 35, 40, math.inf]}) as study:
        rows = study.evaluate_batch(first, second, seeds=seeds, threads=4)
        # Optional NumPy route; contiguous inputs avoid conversion copies.
        array = study.evaluate_array(first, second, seeds=seeds, threads=4)
        labels = study.output_labels
        resolved = study.resolved_configuration()
```

Inputs `first` and `second` have four components `(E,px,py,pz)` per point,
in GeV. Changing cuts or bins builds a new measurement without repeating
evolution. Changing physical or numerical evolution settings requires a new
calculator. Histogram edges still partition the full fiducial support, including
overflow where applicable.

## Production tensors and reusable responses

```cpp
std::vector<DecayResponseNode> nodes{{{100,0,0,0},42}};
DecayAngularGrid angular{16,32};
auto table = DecayResponseTable::build(decay, measurement, nodes, angular, 4);
table.save("responses.lfr");
auto cached = DecayResponseTable::load("responses.lfr", decay, measurement, angular);
auto prediction = cached.response(0).contract(production_tensor);
```

`production_tensor` must follow the documented Cartesian frame, normalization
and supported spin conventions at `cached.node(0).parent`. The table has explicit
nodes and introduces no interpolation. New production variations reuse the
stored tensors. Different cuts, binning, scales, thresholds, sampling controls,
angular nodes, library version or numeric precision invalidate loading.

The example under `examples/response_table/` builds and contracts a small table.
Angular quadrature and hard-sampling accuracy remain independent of persistence.
The native DYTurbo client ({download}`integration guide <../dyturbo_integration.md>`) now runs QCD production
through NNLO+NNLL. The differential NLO+NLL W client also exports reusable QCD
records, and the inclusive Z mass adapter uses a cumulative two-leg loss table.
A validated interpolating spin-response grid and fast correlated finite-term
integration remain separate work; no NNLO+NNLL plus QED phenomenology result is implied.

## Reproducing timing measurements

With examples enabled, run `build/examples/leptonff_decay_timing 256`. It emits
CSV for W/Z photonic and leptonic NLL, original and finer bins, construction,
compiled point actions, four-worker batches, angular responses and cached tensor
contractions. The input point, seeds and numerical settings are in
`examples/decay_timing.cpp`. Response timings use 8x16 angular nodes solely as a
performance probe. Results depend on hardware, occupancy and phase space and
are not numerical-accuracy tests. No large QCD run is launched by this example.

`examples/decay_timing_legacy.cpp` uses only the 0.26 public API, so the same
source can benchmark both releases. It is built as `leptonff_decay_timing_legacy`.
Use its `events` and `leptonic` modes for the recorded comparison; `stats`,
`hard` and profile modes isolate individual paths. The initial measured results
and their settings are in `benchmarks/architecture_0_27.json`.
