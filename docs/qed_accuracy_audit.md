# Accuracy contract for the W/Z comparison

Version 0.24 adds an optional [unitary-conditioned photonic NLL
mode](photonic_nll.md), separate from the photonic LL curve in the frozen
v0.23 campaign described below. The fixed-alpha NLL mode includes the
photonic time-like NLO kernel; pairs and running remain deliberately absent.
The production plots below have not been regenerated with this new curve.

The target is NLO QCD production combined with perturbative final-state QED
fragmentation and a finite, spin-correlated order-alpha decay match. Strong
corrections to production and strong evolution of fragmentation functions are
different accuracy choices. Mixed QCD–QED splitting kernels and a fitted
nonperturbative quark fragmentation boundary are not prerequisites to the
pure-QED NLL benchmark. Earlier status entries that presented them as such
were too broad.

## Relation to the Frixione implementation

[Bertone et al., arXiv:2207.03265](https://arxiv.org/pdf/2207.03265),
sections 3–5, evolve leptons, quarks and photons with pure-QED kernels through
order alpha squared; the gluon PDF is omitted. They match across fermion
thresholds and discuss hadronic vacuum polarisation in the electromagnetic
input. That calculation does not require a fitted hadronic fragmentation
function. It supplies an ISR/space-like reference, rather than the final-state
evolution used here. In particular its zero radiatively generated PDF at a
mass threshold must not be substituted for a time-like FF threshold map.

[Frixione, arXiv:1909.03886](https://arxiv.org/abs/1909.03886), section 6,
derives perturbative initial conditions for fragmentation into electrons and
photons. This is the relevant distinction from fragmentation into hadrons.
The lepton mass sets the perturbative collinear boundary.

## Audited implementation

| Ingredient | Implementation and evidence | Accuracy limit |
| --- | --- | --- |
| Time-like kernels | `src/nlo_lepton.cpp`, `src/nlo_qed.cpp`; direct QED lepton formulas and independently checked time-like continuation for photon channels | Retains a P0 + a² P1, with a=alpha_MSbar/(2 pi); no mixed strong kernels implied |
| Flavour factors | `src/flavours.cpp`; produced-colour factors separated from closed-loop S2/S4; arbitrary-charge momentum/valence checks | Leptonic preset contains e, mu, tau; a separate quark prescription must be named explicitly |
| Boundary | Perturbative one-loop MSbar d1, including plus and delta terms | Unknown two-loop finite boundaries are outside NLL |
| Coupling | Two-loop pure-QED running with active-flavour charge sums; leptonic OS/MSbar conversion at the starting scale | An OS input must not be reused numerically as an MSbar input without conversion |
| Thresholds | Finite time-like heavy-row and photon maps; independent Mellin checks | One-loop maps; finite heavy-row term survives at the mass threshold |
| Hard decay | Muon and photon hard-parent channels, charged-W emission, parity-even neutral-current counterpart | Unitarily normalised, factorised decay, leading power in the muon mass; not full NLO electroweak production |
| Matching | Actual lattice C1 subtraction and continuum restoration; common alpha at original hard mass; consistent retained C2 | Cubic remainders check re-expansion; not a claim to complete order-alpha-squared constants |
| Measurements | Observed radiated momenta, exact endpoint ownership, conservative fraction cells, signed host-event covariance | Inclusive primary-charge W muons and distinct-parent opposite-charge Z pairs |
| Production | MCFM NLO QCD matrix elements, signed real/dipole grouping and common radiation streams | Generates factorised mixed terms, not the complete fixed-order mixed QCD–EW correction |

The source-level kernel and threshold derivations are documented in
[the lepton audit](nlo_kernel_audit.md), [photon audit](nlo_photon_kernels.md),
[flavour contract](flavour_threshold_matching.md) and
[hard matching contract](coupled_lab_matching.md). In particular, the
space-like multi-flavour reference supplies charge/colour bookkeeping; it is
not substituted for the time-like P1 expressions.

The targeted kernel, flavour evolution, independent flavour reference,
independent NLO kernel reference, coupled laboratory matching, fraction-cell,
spectrum, covariance and paired-comparison tests pass (11/11). The log is
`results/qed_accuracy_audit_v023.log`. The last full v0.22 suite passed 60/60;
this comparison changes campaign/analysis/documentation code and adapter
control routing, not the core physics library. Internal consistency and source checks do not establish a
numerical precision claim for the final spectra.

## Three separately labelled comparisons

1. **Photonic:** fixed alpha, no pairs in either tool. The analytic result
   retains LL evolution and at most one finite insertion globally. This is
   the direct coupling/pair-aligned PHOTOS comparison; it is not full NLL.
2. **Leptonic NLL:** e/mu/tau, running MSbar coupling, pure-QED P0/P1 and finite
   thresholds. This tests additional running/pair/NLL content against the
   photon-only reference. A PHOTOS difference is diagnostic because the
   contents differ.
3. **QED + HVP, pointlike prescription:** the existing v0.22 quark/HVP model
   with `--qcd-off`, zero light-quark boundary at the declared infrared
   matching scale, finite heavy thresholds and the sourced fermionic MSbar
   alpha anchor. The gluon is an inert spectator. Strong FF evolution is
   disabled; the sourced alpha_s-dependent high-scale input conversion is
   retained. This is a declared pure-QED evolution prescription with a
   hadronic boundary assumption, not a reproduction of eMELA's mass choices
   or a fitted physical hadronic fragmentation prediction. Its matching-scale
   and HVP-input sensitivity can be measured separately. HVP input alone does
   not fix real hadronic-pair fragmentation.

All use the same NLO QCD production setup and bare-muon cuts. PHOTOS 3.56
has exponentiation and interference enabled, the W matrix-element correction
enabled, and the Z correction disabled. Its massive Born map differs slightly
from the analytic massless map. Thus the Z difference also probes different
finite radiation approximations. The export includes
`(PHOTOS post - PHOTOS Born) - (analytic matched - analytic Born)` so that
the Born-map difference can be identified separately.

## Numerical and phenomenology requirements

The comparison campaign increases calls tenfold relative to the earlier
0.05-scale controls. Independent seeds, doubled grids, wider log-loss domains
and tenfold smaller geometric cutoffs are retained as separate runs. The HVP
domain changes from 64 to the coupled backend's supported maximum of 80.
Fixed-resolution tail controls double both span and bins for the
photonic/leptonic branches; the HVP branch compares span 32 with 64 at the
same resolution. Both HVP domains satisfy the weighted wrap/tail bounds.
The seed samples can be averaged with their independent event covariances.
For the W+ and Z photonic central results, angular-16 estimates replace
angular-4 for independent seeds 1045 and 2089. W− uses ten independent
angular-16 seeds at QCD dipole alpha=0.01 in all three modes, following a
production-variance study. The original alpha=1 runs remain numerical
controls. The ensemble guard checks identical physics and independent seeds;
it permits the technical subtraction-parameter change. Estimates are averaged
with equal weights; correlated original counterparts are not added as extra
independent samples. Independent-seed covariance is retained as a diagnostic.
Refinement shifts compare to the nominal seed, not to the two-seed average.
No cross-run event products are available: separate-run quadrature errors and
paired ten-iteration diagnostics are reported separately. A small noisy shift
does not certify a sub-permille numerical limit.

The requested plots are W± pT and folded |eta|, and Z/gamma* pair mass.
Every plot uses the observed muon momenta. Infinite pT overflow bins remain
in rates, normalisation, covariance and CSV exports; no finite density is
invented for them. Absolute PHOTOS/analytic ratios and normalised-shape
differences retain the full within-run cross-channel covariance. Only the
photonic ratio has aligned coupling/pair content; NLL ratios are diagnostic.
Ratios are
withheld if their signed denominator is compatible with zero at three sigma;
the underlying bin integrals are retained.

The benchmark is pp at 13 TeV, pT(mu)>25 GeV, |eta(mu)|<2.4,
60<m(mu mu)<120 GeV for Z, with no W missing-pT or mT cut. It uses the
already installed central NNPDF31 NNLO set on both NLO QCD calculations;
this preserves a common input for the tool comparison, not PDF-order or
PDF-uncertainty certification. Precision data work still requires a selected
measurement, matching prompt/extra-pair and fiducial definitions, production
PDF/alpha_s uncertainties, and recoil resummation where fixed-order NLO QCD
is insufficient near Jacobian or low-qT regions. None requires us to first
implement the optional combined QCD–QED fragmentation programme.

The domain audit found and fixed a bridge bug: the photonic branch previously
ignored the requested span and always used 40. The corrected bridge honours
the request and logs the resolved grid. The runner records it and checks
agreement; the analyser rejects photonic domain controls without this
evidence. Default old/new W±/Z host outputs agree bitwise in 48 checks,
including all rates, spectra, raw iterations and event products
(`results/qed_adapter_regression_v023/report.json`). The core library, host,
PHOTOS and Fortran measurement/covariance sources are unchanged. The actual
campaign retains earlier default runs and uses corrected domain controls
from `results/qed_domain_v023`, with explicit equivalence evidence.

Fresh reproduction (use a new output directory):

```sh
python3 adapters/mcfm/run_comparison.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-qed-v023/mcfm_leptonff \
  --results results/my_qed_comparison --models photonic leptons qed_hvp \
  --cases nominal seed grid domain tail cutoff --call-scale .5 --jobs 3
python3 adapters/mcfm/run_comparison.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-qed-v023/mcfm_leptonff \
  --results results/my_qed_comparison --models photonic \
  --cases angular angular_seed photos_cutoff --call-scale .5 --jobs 3
python3 adapters/mcfm/run_variance_ensemble.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-qed-v023/mcfm_leptonff \
  --results results/my_variance_ensemble --call-scale .5 --jobs 6
python3 adapters/mcfm/analyse_comparison.py results/my_qed_comparison \
  --cases nominal seed grid domain tail cutoff \
  --extra-photonic-cases angular angular_seed photos_cutoff \
  --photonic-central angular --photonic-seed angular_seed \
  --central-runs results/my_variance_ensemble/central_runs.json \
  --output results/my_qed_comparison/report.json
python3 adapters/mcfm/plot_comparison.py results/my_qed_comparison/report.json \
  --output results/my_qed_comparison/plots
```

The [completed benchmark](../benchmarks/qed_comparison_reference.md) records
144 passing main checks, achieved bin errors and the remaining grid/geometric
precision limits. Passing internal checks does not close those precision or
measurement-definition requirements.
