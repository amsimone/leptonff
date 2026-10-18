# Three introductory calculations

## 1. Fragmentation moments and bins

The complete source is `examples/fragmentation/main.cpp`. Build its standalone
project against an installation:

```sh
cmake -S examples/fragmentation -B build-example -DCMAKE_PREFIX_PATH="$PWD/install"
cmake --build build-example
./build-example/leptonff_fragmentation_example
```

The essential fixed-order calculation is:

```cpp
#include <leptonff/core.hpp>
using namespace leptonff;
const auto ff = fragmentation_order1(
    Species::muon_minus, Species::muon_minus, QEDParameters{});
const Real momentum = ff.moment(2);
const Real endpoint_bin = ff.integral(.9L, 1);
```

The default one-family MSbar input has $\mu_0=m_\mu=0.1056583755$ GeV,
$\mu_F=91.1876$ GeV and $\alpha=1/137.035999084$.
The N=2 result is approximately `0.976868`; the precise value and tolerance
are distributed in `benchmarks/canonical/order1.json`. This fixed-order part
takes well below a second on an ordinary laptop; the example also evaluates
a continuum coupled NLL moment with explicit leptonic on-shell conversion.

Bin conventions are `[lower,upper)`, with the endpoint included in the final
bin when `upper=1`. This matters for delta contributions. Full number integrals
of soft-divergent channels, such as muon-to-photon, are not finite.

## 2. Convolve with a coefficient function

`examples/convolution/main.cpp` evaluates

$$
(D\otimes c)(x)=\int_x^1\frac{dz}{z}D(z)c(x/z),\qquad c(x)=x(1-x).
$$

```cpp
const Function coefficient = [](Real x) { return x*(1-x); };
const Real value = ff.convolve(coefficient, .5L);
```

The convolution includes the delta and whole-plus pieces. Its domain is
$0<x<1$, with a coefficient smooth on $[x,1]$. Known discontinuities require
the lower-level `Distribution::apply` interface with explicit breakpoints.
This example takes well below a second and has no generator dependency.

The same calculation in Python is:

```python
from leptonff import Fragmentation
with Fragmentation() as ff:
    print(ff.moment(2))
    print(ff.integral(.9, 1))
    print(ff.convolve(lambda x: x*(1-x), .5))
```

## 3. Matched W/Z response

`examples/cpp/main.cpp` creates an immutable `DecayCalculator`, evaluates one
W Born point and contracts a response with a production spin tensor. Build it
using its standalone CMake project, or run
`build/examples/leptonff_matched_example` in the default example build.

The example uses 128 lattice bins, four hard samples and a 4-by-8 angular
quadrature to keep it small. These are demonstration controls, rather than a
phenomenological precision configuration. It prints Born acceptance `1` for
its supplied W point, the matched response and a tensor contraction. Runtime
is typically seconds, depending on compiler and machine.

For a simple Python W point:

```python
from leptonff import Decay, DecayConfig
config = DecayConfig(grid_bins=128, hard_samples=4)
with Decay(config) as decay:
    prediction = decay.evaluate([50,40,0,30], [50,-40,0,-30], seed=42)
    print(prediction.total.born, prediction.total.matched)
```

Supported observable names are W: `pt`, `abseta`, `mt`; Z: `mass`, `minpt`,
`maxeta`. Bins must cover their complete configured fiducial support. `mass`
denotes $m_{\ell\ell}$, `minpt` the smaller lepton transverse momentum, and
`maxeta` the larger absolute pseudorapidity.

For scans, prepare `QEDXSpaceEvolution` or `DecayCalculator` once and reuse it.
The current implementation supports in-memory reuse. Portable serialization
and interpolated W/Z response tables remain future extensions.
