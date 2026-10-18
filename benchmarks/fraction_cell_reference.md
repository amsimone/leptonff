# Conservative fraction-cell benchmark

Version 0.19 improves the numerical measurement of the coupled one-family
W±/Z hard match. Resolved weights are integrated over conservative signed
cells; the endpoint atom remains exact. Evolution, hard physics and input
conventions are those of the [v0.18 benchmark](coupled_lab_reference.md).
See the [cell contract](../docs/conservative_fraction_cells.md) for the exact
number/energy conservation and hyperbolic pair-cut integrals.

## Run contract

The nominal grid has 4096 bins over log-loss span 40, minimum fraction 1e-4,
and four hard samples per QCD entry. The 13 TeV MCFM NLO QCD, bare-muon cuts,
one-family MSbar input and diagnostic PHOTOS setup are unchanged. New runs
have the explicit `--fraction-cells` flag. Atomic controls use the same new
executable with that flag inactive.

Each process has nine runs: nominal, atomic, zero coupling, 8192/16384 bins,
QCD dipole alpha=.1, geometric exclusions 1e-6, and cell/atomic runs with
independent seed 2089. The ordinary seed is 1045. Every comparison retains
signed QCD real/dipole grouping and common radiation streams. The optional
raw iteration audit includes all 33 rates.

Statistics use common equal-iteration estimates with directly accumulated
higher-order differences. Per-histogram native inverse-variance averages are
not subtracted to construct an error. Paired grid diagnostics use ten chunks
per QCD contribution, checking identical Born controls at every iteration;
they have no event-level cross-products or rare-tail guarantee.

## What is tested

The 53-test core suite checks exact signed cell overlaps and reciprocal
prefixes against literal integration, number/energy conservation, both
strict endpoint atoms, mixed atomic/cell pair integrals and symmetry, and
resolved CDF continuity. An independent D(z)=1 fixture verifies second-order
CDF convergence to its known continuum integral. This smooth fixture is not
a general error theorem for physical endpoint singularities.

W±/Z C1/C2 tests use the actual cell measure in the evolved term and overlap,
leaving cubic remainders after subtracting through second order. Zero
coupling reduces exactly to the Born measurement. The GCC host compiler
also passes the cell and coupled tests. Active-cell PHOTOS bridge checks
cover all 33 outputs at 1000 events per process, including exact common-seed
repeatability and conservation/radiation checks.

Host checks include inactive-mode bitwise reproduction of every v0.18 rate,
unchanged continuum NLO and PHOTOS columns, zero identities, grid/source
invariants, QCD subtraction and geometric-cutoff controls, and independent
seed consistency. Measured finite-grid cell/atom differences are diagnostics,
not expected exact zeros: the two numerical measures differ beyond O(alpha).
The full report is `results/fraction_cells_host.json`.

## Deterministic grid response

`adapters/mcfm/check_cell_response.cpp` evaluates boosted fixed W and Z decay
configurations with an unpolarised density and two fixed hard-radiation
points, at 4096–32768 bins. It asserts the requested fraction mode is active.
This removes QCD sampling noise from the grid comparison. The reported
quantity is a **dimensionless fixed-point higher-order integrand**, not an
integrated decay acceptance or a hadronic cross section. The raw CSV,
component closure and source/library/executable hashes are in
`results/cell_fixed_response/`.

For the first hard point, values below are multiplied by 10^6:

| Bins | W atomic | W cells | Z atomic | Z cells |
| ---: | ---: | ---: | ---: | ---: |
| 4096 | -55.7870 | -55.4279 | -146.3092 | -146.8930 |
| 8192 | -55.4308 | -55.3534 | -148.5968 | -146.4569 |
| 16384 | -55.2698 | -55.3329 | -147.1606 | -146.3505 |
| 32768 | -55.3345 | -55.3274 | -146.2904 | -146.3216 |

The selected cell responses show a smoother trend. At the last doubling,
the changes are about twelve times smaller for W and thirty times smaller
for Z than the atomic changes. The second hard point has the same qualitative
behaviour. These checks do not certify a uniform phase-space or total-rate
error bound.

## Reproduction

```sh
python3 adapters/mcfm/build_adapter.py /path/to/MCFM-10.3.1/build-gcc \
  --core-build build-host --output build-mcfm-cells \
  --photos-build /path/to/photos-build --iteration-audit
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-cells/mcfm_leptonff --threads 1 \
  --output results/my_wplus_cells \
  --coupled-hard --fraction-cells --hard-samples 4
python3 adapters/mcfm/run_cell_controls.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-cells/mcfm_leptonff --results results/new_cells \
  --cases nominal atomic zero grid8192 grid16384 alpha01 slivers seed2089 atomic2089
```

Repeat the controls for W- and Z. Every run directory must be new. Audit a
new result root against the retained v0.18 reference with:

```sh
python3 adapters/mcfm/analyse_cells.py --results results/new_cells \
  --reference-results results --output results/new_cells_audit.json
```

The first executable source is frozen in `results/cell_source_v019/`, core
digest `bd5e4022975774d20b33d91d5ede02a9a23f9626a33b2a32b4b73a3b0d94db4c`.
The final extra smooth-continuum regression test is frozen in
`results/cell_source_v019_final/`, core digest
`0bb1eefc07874046752d35706e896d6d5e508dc58661fe2ce9b2dfdbec8de8bd`.
Only that test differs; the production physics library is unchanged between
those snapshots. Run manifests distinguish executable and working-tree
source hashes. Final analysis/diagnostic scripts are archived separately.

## Remaining precision work

The observed variance/CPU gains vary by seed and observable. The exact
endpoint retains a physical Born cut, and rare QCD/angular contributions
remain. Neither this reconstruction nor a small fragmentation-scale band
provides a complete numerical/statistical budget. Physical flavours and
thresholds, coupling input conversion, hard differential spectra and an
experimental pair-selection prescription remain separate work.

## Completed host audit

All **36 checks pass** across 27 completed host runs, with no diagnostic
flags or integration warnings. The maximum required subtraction/cutoff pull
is 1.57 combined standard deviations. Inactive cell mode reproduces every
v0.18 rate bitwise, including PHOTOS and the coupled higher-order columns.
The active mode retains exactly the same continuum NLO and PHOTOS columns.

The final core suite passes **53/53** in 152.69 seconds, in
`results/cell_final_ctest.log`. Active-cell bridge checks and both GCC core
checks pass; their reports are `results/cell_smoke/report.json`,
`results/cell_gcc.log` and `results/cell_gcc_final.log`.

For the directly accumulated matched-minus-continuum-order1 rate, observed
variance gains (atomic variance divided by cell variance) and CPU efficiency
gains are:

| Process | Variance gain, seed 1045 | Variance gain, seed 2089 | CPU efficiency gain, seed 1045 | CPU efficiency gain, seed 2089 |
| --- | ---: | ---: | ---: | ---: |
| W+ | 2.523 | 1.091 | 2.483 | 1.084 |
| W- | 1.087 | 1.301 | 1.066 | 1.286 |
| Z/gamma* | 1.184 | 1.204 | 1.233 | 1.153 |

CPU efficiency is `(variance_atomic*CPU_atomic)/(variance_cells*CPU_cells)`;
all runs use one thread. These are sample-specific measurements, not
uncertainties on the gains or guarantees over phase space. In particular the
larger W+ gain of the first seed is not universal.

Nominal cell-measured common-iteration rates are:

| Process | Matched [pb] | Matched minus continuum order1 [pb] |
| --- | ---: | ---: |
| W+ | 5383.2579 ± 22.8814 | -0.08337 ± 0.77130 |
| W- | 4126.3895 ± 36.4840 | +0.54979 ± 0.58664 |
| Z/gamma* | 725.2350 ± 1.8129 | +0.02982 ± 0.08049 |

Grid changes below are relative to 4096 bins. Their errors are limited paired
iteration-chunk diagnostics. They are not certified numerical error bounds.

| Process | Atomic 8192 shift [pb] | Cell 8192 shift [pb] | Atomic 16384 shift [pb] | Cell 16384 shift [pb] |
| --- | ---: | ---: | ---: | ---: |
| W+ | +0.2906 ± 0.7605 | +0.2218 ± 0.2257 | +1.1164 ± 1.0936 | +0.4624 ± 0.4176 |
| W- | +0.1552 ± 0.2091 | +0.0480 ± 0.1276 | +1.1633 ± 0.8833 | +0.5920 ± 0.5167 |
| Z/gamma* | -0.0716 ± 0.0464 | -0.0086 ± 0.0178 | -0.1181 ± 0.0548 | -0.0291 ± 0.0322 |

The grid comparisons are less noisy and their measured changes are smaller,
but the residuals remain comparable to the small higher-order terms. A full
continuum-rate error and rare-tail control have not been established. The
finite-cell representation is therefore an explicit opt-in numerical method,
with the atomic reference retained for further convergence studies.

Final runner, audit and fixed-response sources, plus report/library/binary
hashes, are archived in `results/cell_analysis_v019/`. No earlier scientific
result is relabelled with the final regression-test source digest.
