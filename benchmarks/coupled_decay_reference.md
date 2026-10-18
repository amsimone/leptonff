# Coupled fiducial hard-matching reference

Version 0.11.0, Release, AppleClang 21, 53-bit Real. All **40 CTest checks pass**.
Pure unpolarised vector at rest, one muon family, inclusive identified
opposite-charge pairs from distinct hard parents. No extra-pair veto.
Normalisation is the massless inclusive vector Born rate, not a physical DY
cross section. The [contract](../docs/coupled_fiducial_matching.md) specifies
the matching, signed quadrature and approximation.

Inputs: M=91.1876 GeV, m=mu_0=0.1056583755 GeV, explicit
alpha_MSbar(mu_0)=0.0072973525692838015, pT>25 GeV, |eta|<2.4,
60<m_pair<120 GeV. Two-loop running supplies alpha_MSbar(M)=0.007374692627489.
The input is a model MSbar coupling, without an on-shell conversion.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure --parallel 4
./build/leptonff_coupled configs/reference.conf 8192 8 --scan \
  > results/coupled_decay_fine_scan.json
python3 python/check_coupled_decay_reference.py build/leptonff_coupled \
  configs/reference.conf --input results/coupled_decay_fine_scan.json \
  --output results/coupled_decay_reference.json
```

The saved checked scan includes source/executable/configuration/input hashes
and an independently checked 256-bin/two-panel fixture. `--dump` exports the
hard prefixes, FF/coefficient arrays, interpolation fractions and geometry
matrices. A separate initial scan used 4096 bins/four panels.

## Numerical convergence

| Bins | Log span | FF panels | Hard controls | Matched rate |
| --- | --- | --- | --- | --- |
| 4096 | 40 | 4 | baseline | 0.806898845499 |
| 4096 | 40 | 8 | baseline | 0.806898937493 |
| 8192 | 40 | 8 | baseline | 0.806900029940 |
| 8192 | 40 | 16 | baseline | 0.806900020737 |
| 16384 | 40 | 16 | baseline | 0.806900274387 |
| 32768 | 80 | 16 | baseline | 0.806900274387 |
| 32768 | 80 | 16 | refined | 0.806900275699 |

Baseline geometry controls: soft and collinear slivers 1e-6; absolute/relative
nested GK tolerances 2e-6. Refined controls divide both slivers and tolerances
by four. Hard source projection tolerances are 1e-9. The full O(alpha) overlap
is restored at each point, with the same geometry controls as its continuum
target.

Final scan shifts:

| Refinement | Absolute rate shift |
| --- | --- |
| FF panels 8 -> 16 at 8192 bins | -9.2035e-9 |
| Grid 8192 -> 16384 | +2.5365e-7 |
| Span 40 -> 80 at fixed finest spacing | +6.66e-16 |
| Hard slivers/tolerances divided by four | +1.3121e-9 |

The last grid shift dominates these numerical comparisons. They support
roughly a few 1e-7 to 1e-6 absolute numerical control for this fixture, rather
than every printed digit or a universal bound for sharp/narrow cuts.

At the final point, weighted FF wrap/tail bounds are 1.04e-26 and 7.49e-36;
the weighted step-doubling shift bound is 5.44e-9. These are energy-weighted
solver diagnostics, not total fiducial error bounds. The hard projections use
causal sums and have no added periodic-wrap approximation. Signed cubic FF
compression preserves truncated moments through degree three to roundoff;
its sharp-cut interpolation error is measured by the panel refinement.

## Components and scale variation

At the final central point:

| Component | Contribution to rate |
| --- | --- |
| Born | 0.836268494675 |
| Raw coupled Born-FF term | 0.796739138289 |
| Muon-parent hard term, including a(M) | +0.010168162661 |
| Photon-parent hard term, including a(M), both primary charges | -0.000006829268 |
| Full O(alpha) overlap restoration | -0.000000195984 |
| Matched total | **0.806900275699** |
| Continuum full O(alpha), expanded in a(mu_0) | 0.806597111540 |

The raw hard coefficients before a(M) are 7.522160351756 (muon projection),
1.141041419698 (muon geometry), -0.007773564116 (photon projection) and
0.001955076746 (photon geometry). Photon entries already sum the two primary
charges. The explicit second coefficient of the numerical prescription is
222.485358145 in a(mu_0)²; constant two-loop boundary/H2 terms are absent.

| mu_frag/M | Matched rate |
| --- | --- |
| 0.5 | 0.806915004535 |
| 1 | 0.806900275699 |
| 2 | 0.806869340189 |

The hard coupling stays at M. The retained alpha² mass logarithm cancels in
the scale shift; unsupplied NNLL terms leave a nonzero residual band. This is
not a demonstrated physics-uncertainty band for physical DY, and the old
fixed-alpha scalar LL number is not an aligned LL comparator.

## Independent checks

The coarse fixture is validated with source cells in z, nonperiodic causal
power-series evolution, analytic hard moments/endpoints, literal pair sums,
fixed cubic interpolation and direct allowed-azimuth geometry. Production
instead uses GK cells, coupled FFT/RK evolution, direct endpoint hat actions,
prefix-sum Born measurements and spherical-cap intersections.

| Check | Maximum absolute residual |
| --- | --- |
| Truncated evolved FF weights | 1.52e-13 |
| First/second FF coefficient arrays | 1.01e-8 |
| Hard source projection | 5.08e-9 |
| Independent finite-geometry contractions | 4.23e-6 |
| Geometry 96 -> 192 quadrature refinement | 1.61e-5 |
| Full matched fiducial rate | **5.13e-9** |
| Actual first overlap coefficient | 3.60e-7 |
| Second coefficient with independent sharp-cut geometry | 3.93e-3 |
| Second coefficient with the separately audited production geometry | 3.73e-8 |

The independent second coefficient is limited by angular quadrature; its
96-to-192 refinement is 1.59e-2, equivalent to about 2.2e-8 in the rate. Keeping
the audited geometry fixed isolates the much tighter expansion-algebra check.
The final rate agreement benefits from correlated overlap subtraction and is
not a claim of 1e-9 precision in every ingredient.

C++ tests separately check the cubic small-coupling remainder after removing
the explicit first and second coefficients, O(alpha) scale cancellation,
absence of a retained alpha² mass logarithm in scale variation, common hard
coupling with the moment reference, zero coupling and invalid configurations.
The existing massive-current/fiducial checks establish the leading-power
fixed-order target. No external full NLL generator comparison has been made.

The restricted one-family fiducial rate is now joined and validated. Binned
matched output, physical spin/current/flavour/threshold content, a full
uncertainty budget and an NNLO production adapter remain outstanding.
