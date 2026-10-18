# LL-dressed hard-decay reference

Version 0.6.0, AppleClang 21, 53-bit Real. Both existing configurations work:

```sh
./build/leptonff_dressed configs/reference.conf
./build/leptonff_dressed configs/vector_mass.conf
python3 python/check_hard_dressed_reference.py build/leptonff_dressed configs/reference.conf \
  --output results/hard_dressed_reference.json
```

The production angular calculation is slower than the additive driver because
it integrates the hard remainder for pairs of FF quadrature nodes. Optional
controls are `config bins LL_cutoff log_span FF_panels`. The default uses 16
panels, then 32; each positive panel has up to two Gaussian nodes, in addition
to the exact endpoint atom and a node retaining signed FFT roundoff. The driver
also doubles the lattice, halves its soft cutoff, then refines the hard-geometry
controls. Work in the angular remainder grows quadratically with node count.

For the pT>25 GeV, |eta|<2.4, 60–120 GeV fixture, relative to massless vector Born:

| Stage | Matched rate |
| --- | ---: |
| Full O(alpha) | 0.806597111174 |
| Undressed hard + LL/boundary, previous driver | 0.807232541918 |
| Dressed baseline | 0.806839089814 |
| Doubled FF quadrature | 0.806836330287 |
| Refined LL grid | 0.806836334749 |
| Refined LL soft cutoff | 0.806836607035 |
| Refined hard geometry | 0.806836606334 |

The final hard coefficient is 8.687532127433. It equals the undressed value
9.028441603971 plus a projected dressing shift -0.257479099853 and a finite
geometry shift -0.083430376685. These coefficients multiply alpha/(2*pi).

| Separate numerical change | Rate shift |
| --- | ---: |
| FF quadrature refinement | -2.760e-6 |
| LL grid refinement | +4.463e-9 |
| LL soft-cutoff refinement | +2.723e-7 |
| Hard soft/angular/tolerance refinement | -7.004e-10 |

The displayed digits permit reproduction; they are not a twelve-digit accuracy
claim. FF quadrature currently dominates these numerical shifts. Convergence
need not be monotonic; refine further for a tighter target or different cuts.
These shifts are numerical checks, not estimates of missing physical orders.

For the angular-inclusive 60–120 GeV mass window, the geometry term vanishes.
The final dressed rate is 0.985766562232, versus 0.985462761542 with an undressed
hard term and 0.986414026151 at fixed order. Its grid and LL-cutoff shifts are
+2.140e-8 and -1.465e-7; FF quadrature has no effect.

The independent coarse Poisson check verifies signed hard cells and their
convolution. Direct azimuthal 32/64-point angular rules check the dressed
geometry; the fine coefficient differs from production on that coarse fixture
by 8.14e-6. Polynomial moments through z^3 agree with the independent lattice
weights. The mass-only check additionally reduces the two muons to U(2*t).
The C++ coupling scan establishes that dressing changes the old additive
result first at alpha squared; zero evolution recovers full O(alpha).

Production manifests record all controls, endpoint/node weights, hard cell
heads, FFT diagnostics, wrap/tail bounds and source/executable/configuration
hashes. To validate an already saved driver output without repeating the costly
angular run, pass `--input path/to/raw.json` to the Python checker. The input
file's SHA-256 is recorded separately.

All 21 CTest checks pass. This remains fixed-alpha, one-family, unpolarised
vector decay at rest. The NLO splitting kernel, running/mixing and physical
spin/boost extensions are absent. [Derivation and scope](../docs/hard_dressing.md).
