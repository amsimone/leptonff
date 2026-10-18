# W/Z photonic hard-dressing benchmark

Version 0.17. This implements the next finite-dressing increment on the NLO
QCD host: one MSbar boundary or spin-correlated H1 insertion, dressed with LL
on all charged hard-event legs. See the [operator contract](../docs/photonic_lab_dressing.md).
It is **not full NLL**: P1/running and physical flavour/pair content are absent
from this comparison branch.

## Inputs and rates

The v0.16 inputs are retained: MCFM 10.3.1 NLO QCD, 13 TeV pp,
NNPDF31_nnlo_as_0118 member 0 as the available validation input, fixed MW/MZ
QCD and fragmentation scales, alpha=1/137.035999084 and m_mu=0.1056583755 GeV.
W is on shell. Z/gamma* generation includes parents above the 120 GeV observed
mass boundary. Both observed muons require pT>25 GeV and |eta|<2.4, with
60<m_mumu<120 GeV for Z and no neutrino-pT or mT cut for W. No recombination.

The LL grid has 4096 bins over log-loss span 40, with a half-cell unresolved
cutoff. Four geometric samples are used per QCD entry; real and mapped
counterevents share the same original-real radiation seed. Finite geometric
soft and 1-|c| exclusions are 1e-5. PHOTOS 3.56 is exponentiated, at physical
alpha with cutoff 1e-5, interference on, W ME on and Z ME off.

All estimates below use **common equal-iteration weights** within each QCD
contribution, including zero iterations. Errors combine the stored raw
within-iteration variances; the ten-chunk covariance is also retained.
They are integration diagnostics, not a certified rare-event or physics
uncertainty budget. The dressing difference is accumulated directly.

| Process | Dressed matched rate [pb] | Paired dressing shift [pb] | Paired PHOTOS minus dressed [pb] |
| --- | ---: | ---: | ---: |
| W+ | 5382.98 ± 22.92 | -1.320 ± 0.726 | -2.46 ± 6.40 |
| W- | 4125.86 ± 36.60 | -0.469 ± 0.590 | -37.43 ± 32.33 |
| Z/gamma* | 725.145 ± 1.817 | -0.1704 ± 0.0637 | -2.058 ± 0.987 |

The large W- common-weight error is retained. It reflects rare-event and
iteration sensitivity already exposed by the previous pilot. Smaller native
per-histogram inverse-variance errors must not be substituted for a
controlled common uncertainty. The roughly 2-sigma Z tool difference includes
the branches' different perturbative content and generator approximations;
it is not an isolated measurement of NLL effects.

The common-weight split of the dressing shift is:

| Process | Boundary [pb] | Projected hard [pb] | Geometric hard [pb] |
| --- | ---: | ---: | ---: |
| W+ | -0.2404 ± 0.0953 | -0.2608 ± 0.1370 | -0.8183 ± 0.7246 |
| W- | -0.2490 ± 0.2789 | -0.0422 ± 0.4033 | -0.1774 ± 0.5767 |
| Z/gamma* | -0.1623 ± 0.0365 | +0.0617 ± 0.0435 | -0.0698 ± 0.0633 |

The direct total error includes correlations and is not the quadrature sum
of the component errors. The raw 26-rate iteration records provide the
corresponding covariance estimate.

## Validation

The core suite passes **51/51 tests**. The new reference checks include:

- Hard source cells against independent kernel quadrature for W and Z.
- Zero-padded hard transport against literal causal convolution.
- Literal one-/two-leg boundary sums, excluding a second boundary insertion.
- Dressed real geometry against separate fragmented spin-event sums for W±/Z.
- Original hard-scale reconstruction for W and an off-shell Z parent migrating
  into the observed mass window; fragmented daughter masses are not used as
  the hard scale.
- Exact identity-evolution/zero-coupling limits, the unchanged NLO overlap,
  alpha-squared hard shifts and explicit boundary C2 extraction.

The optional PHOTOS smoke tests pass with dressing active for all three
processes, including common-seed repeatability and momentum conservation.
Repeating the independent RECOLA 1.5.1 check with the current library retains
all 486 polarised real points, with maximum residuals 3.16e-14 (W-), 3.13e-14
(W+) and 5.59e-15 (parity-even Z). See `results/recola_spin_dressed/report.json`.

Every earlier v0.16 rate column, including continuum NLO and the generator
branch, is reproduced **bitwise** by the new nominal executable. Its new
columns therefore measure only the added finite dressing.

Nominal data are in `results/nlo_{wplus,wminus,z}_dressed_audit/` and
`results/dressed_nominal_iterations.json`. The first `*_dressed/` pilots retain
valid native rate results, but their initial raw-iteration filter omitted new
prefixes; those incomplete raw records are not used for the final audit.
The corrected `*_dressed_audit/` runs reproduce the native rates and retain
all 26 columns. The combined control report is
`results/photonic_dressing_host.json`.

## Reproduction

```sh
python3 adapters/mcfm/build_adapter.py /path/to/MCFM-10.3.1/build-gcc \
  --output build-mcfm-dressed-audit --core-build build-host \
  --photos-build /path/to/photos-build --iteration-audit

python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-dressed-audit/mcfm_leptonff \
  --output results/my_wplus_dressed --threads 1 \
  --photonic-ll --hard-samples 4 --dress-hard
```

Repeat for W- and Z. The controls runner supplies zero coupling, 8192 bins,
QCD dipole alpha=.1 and finite geometric exclusions 1e-6:

```sh
python3 adapters/mcfm/run_dressed_controls.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-dressed-audit/mcfm_leptonff
python3 adapters/mcfm/analyse_dressed.py
```

Each run directory must be new. `analyse_dressed.py` checks native same-stream
invariants and uses the raw common-iteration estimates for subtraction/cutoff
controls and the dressing shifts. Grid shifts are measured changes rather
than residual error bounds.

## Remaining work

Join the spin-correlated hard operator to the coupled NLL evolution with all
required hard-parent channels, pair/charge counting and consistent hard
couplings; extend physical flavours and thresholds; add hard-event differential
spectra; improve the difficult W variance and the joint covariance/error
budget. This increment supplies finite-term dressing for the aligned photonic
tool benchmark and does not close those broader phenomenology requirements.

## Completed controls and source audit

All **21 host checks pass**, with no retained component flags in this new
control set. The required QCD-subtraction and finite-geometric-cutoff
comparisons have maximum displacement 1.41 combined sigma, using common
iteration estimates. The v0.16 statistical flag remains in its historical
report; this does not turn the new limited sample into a precision budget.
Zero-alpha identities are exact for W±/Z. Rate-component closure is within
9.7e-11 pb across every raw production iteration.

Doubling the LL grid from 4096 to 8192 gives these common-iteration shifts:

| Process | Complete dressed rate shift [pb] | Dressing-only shift [pb] |
| --- | ---: | ---: |
| W+ | +0.06007 | +0.04522 |
| W- | +0.14130 | +0.10306 |
| Z/gamma* | -0.06729 | -0.01404 |

The continuum NLO, bare hard and PHOTOS columns are bitwise unchanged by
this grid refinement. Tightening both geometric exclusions by ten passes
all checks; its measured shifts and separate-run errors are retained in the
combined JSON. These are sampling/refinement observations, not certified
residual numerical bounds or complete physics uncertainties.

The final CTest run passes 51/51 in 143.40 seconds
(`results/dressed_final_ctest.log`). Build/run source snapshots are:

- `results/dressed_source_v017/`: first executable/pilot source,
  digest `50c61a4e257df19d2e010a4a98ea3e89f3bcbe46d75fc2bf61003f9b7cf3e058`.
- `results/dressed_source_v017_audit/`: complete-iteration executable source,
  digest `60697daf47f9c2f3df87c4b69cc9944cfeab0f8f1d120d9c37872956cdb8ebfe`.
- `results/dressed_source_v017_final/`: final tested source, with the added
  off-shell hard-scale regression,
  digest `c26ceb11be46f80d1ae4aa8ad59ef25ab1e398126250641f46a1b6510ec12e5d`.

The final additional test does not change the physics library: its library
hash matches both adapter build manifests. Per-run manifests retain the
actual executable, configuration, library and source hashes; no earlier
artifact is relabelled with a later source version.
