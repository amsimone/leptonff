# Validation and numerical precision

## Redistributable canonical benchmark

`benchmarks/canonical/order1.json` records complete physical inputs, an absolute
tolerance and six reference values. They cover the N=1,2,3 same-charge
fragmentation moments, an interior plus-kernel bin, an endpoint bin and a
whole-plus convolution. Values are derived from the analytic formulas stated
in the benchmark's README, independently of the executable under test.

`examples/benchmark.cpp` emits the corresponding numerical values as JSON.
The `canonical_benchmark` test checks them without a production dataset or an
external generator. These are analytic checks, rather than a claim of
phenomenological bin precision.

## Test levels

| Level | Evidence | Purpose |
|---|---|---|
| Analytic | Sum rules, endpoint identities, moments, overlap expansions | Check mathematical invariants and conventions |
| Independent numerical | Python reference implementations and optional external matrix-element/host comparisons | Cross-check implementations over their shared scope |
| Regression | Existing W/Z, threshold, signed-group and replay tests | Preserve established results during changes |
| Numerical convergence | Grid, cutoff, quadrature, hard-sampling and QCD statistics scans | Establish precision for a specified observable |
| Packaging | Relocated C/C++/Fortran consumers and Python/native checks | Verify real installation and ABI behavior |

Core unit checks are enabled with `BUILD_TESTING`. For the complete research
reference suite:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON \
  -DLEPTONFF_BUILD_REFERENCE_PROGRAMS=ON
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --parallel 4
```

Shared-library builds additionally test the Python numerical bindings. Enable
`LEPTONFF_BUILD_FORTRAN` to test an installed compiled module; when only a client
compiler is available, the relocation test exercises the installed module source.

## Establishing phenomenological precision

Quote QCD integration errors separately from fragmentation grid/cutoff errors,
hard-geometry sampling and angular quadrature shifts. Preserve cross-bin and
cross-model covariance when replaying the same signed production groups.
Deterministic quadrature nodes are not independent Monte Carlo events.

Refine each numerical control for the actual cuts and bins. A periodic-wrap
bound alone does not bound sharp-bin discretization error. More bins can require
substantially more production statistics. The package refactor adds no new
high-statistics QCD sample and does not certify a tenfold error reduction.

Historical studies remain under `benchmarks/` and the earlier topical documents
under `docs/`. Some cite large local `results/` files excluded from the source
distribution. They are development evidence, separate from the canonical
redistributable benchmark and the reproducible tests above.

## Complete Z radiator audit

The optional research study in `docs/z_radiator_validation.md` checks the
complete normalized two-lepton mass radiator, including its hard coefficient,
through second order. It compares continuum moments, double-plus cumulative
integrals, mass-log coefficients and coupling expansions with the production
operator. Install NumPy/SciPy before configuring to enable the corresponding
independent CTest; the full study uses finer grids than its CI counterpart.

A separate pinned KKMC adapter supplies fixed-parent-mass comparisons with
cross-order covariance and cutoff/coupling controls. It is not a library or CI
dependency. The study resolves external residuals and keeps them distinct from
the successful implementation checks. Neither those checks nor the QED scale
and matching-prescription variations establish a total phenomenological error.

The later optional study in `docs/mcmule_fixed_order.md` uses pinned McMule
real–virtual/double-real photonic FSR with the inclusive NLO normalization
counterterm. A six-Q scan of strict second-order moments is consistent with
finite constant remainders after subtracting the NLL expansion. Independent
phase-space-volume and FKS/technical controls accompany it. This supports the
mass-log coefficients of the tested inclusive vector observable; sharp endpoint
cumulatives and full finite NNLO matching remain separate precision tasks.
`docs/computation_report.md` inventories physical inputs, evidence and exclusions.
