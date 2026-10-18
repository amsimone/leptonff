# Pure-QED W/Z comparison campaign

The [accuracy audit](../docs/qed_accuracy_audit.md) defines the intended NLO
QCD production plus final-state QED calculation and corrects the earlier
overbroad mixed-kernel/nonperturbative-boundary requirements. The
[completed comparison](../results/qed_comparison_v023/comparison.md) contains
W± pT and folded |eta| and Z/gamma* dilepton-mass figures, absolute ratios,
normalised shapes, bin tables and covariance. Standalone PNG/SVG/PDF exports
are under `results/qed_comparison_v023/plots`.

## Setup and validation

All curves use pp at 13 TeV, pT(mu)>25 GeV, |eta(mu)|<2.4, with
60<m(mu mu)<120 GeV for Z. W has no missing-pT or mT cut. The installed
central NNPDF31 NNLO set is shared by the NLO QCD matrix elements. The direct
PHOTOS comparison uses fixed alpha and no pairs in both tools; the analytic
branch is photonic NLO+LL with one finite insertion globally. PHOTOS 3.56
uses exponentiation/interference, W ME correction on, Z ME correction off.
Leptonic NLL and pointlike QED+HVP curves are separately labelled diagnostics.

- **144/144 main consistency checks pass**, with no main host warning flags.
  These cover 63 baseline/control runs, 30 central W− ensemble runs and shared
  Born/PHOTOS references across the three fragmentation modes.
- **15/15 HVP-input checks pass**, retaining five reduced-call W+ host warnings.
- **48/48 adapter regressions pass:** default W±/Z rates, all 136 spectra,
  raw iterations and event products agree bitwise before/after the domain fix.
- **11/11 focused tests pass**, including independently checked time-like
  kernels/flavour evolution, matching and paired/cross-bin statistics. The
  additional tests cover singular-covariance chi-square diagnostics and
  measuring PHOTOS itself in the PHOTOS-cutoff control.

The 80 C++ source/header files under src/include/apps are unchanged from
the frozen v0.22 implementation. Old and corrected adapters use the identical
core library SHA256
`34d37d38f33d086ae979fbb357a35c1ce9662e2e5891b0655f7cc50bd44c5891`.
The bridge fix honours the requested photonic log span and logs its resolved
value. The analyser rejects domain controls without applied-value evidence.
The runner also enforces the coupled backend's actual span limit of 80.

## Achieved comparison

The photonic central result uses two independent angular-16 seeds for W+
and Z, at ten times the earlier 0.05 validation call budget. A noisy W− bin
was traced to production sampling: its pre-radiation variance was already
large and increasing decay-angle samples did not resolve it. A QCD dipole
alpha study at 1, 0.1 and 0.01 motivated ten independent W− seeds at 0.01,
with angular-16 sampling in all three fragmentation modes. All ten seeds
are retained with equal weights and their independent covariances. The
technical subtraction parameter changes the integration, not the NLO cross
section. Configuration guards check that the physical setup is unchanged.

| Process | Analytic photonic [pb] | Paired PHOTOS − analytic [pb] |
| --- | ---: | ---: |
| W+ | 5383.60 ± 18.34 | −3.83 ± 7.32 |
| W− | 4113.63 ± 8.89 | 4.09 ± 3.22 |
| Z/gamma* | 722.62 ± 1.89 | 0.07 ± 1.48 |

Errors are integration statistics. The paired ratios retain denominator
covariance; these rates/errors are not used as independent quantities when
forming ratios. The problematic W− 25–30 GeV bin now has a ratio error of
0.52%, compared with 5.4% in the original two-seed configuration. Its
independent-seed scatter diagnostic gives 0.73%. W pT bins near the Jacobian
edge still have approximately 1–2% comparison errors. W |eta| errors are
approximately 0.23–0.36%; central ratio departures are at most about 0.5%.

| Normalised photonic/PHOTOS shape comparison | chi-square / covariance rank |
| --- | ---: |
| W+ pT | 5.74 / 5 |
| W+ absolute eta | 3.25 / 3 |
| W− pT | 5.52 / 5 |
| W− absolute eta | 2.41 / 3 |
| Z pair mass | 14.92 / 7 |

The W− independent-seed pT consistency diagnostic is 53.23/54. The Z mass
comparison shows a few-percent difference in some bins; around 90–95 GeV
the absolute ratio departure is −0.18 ± 0.22%. These compare specified
algorithms, including the disabled PHOTOS Z ME correction and different
higher-order approximations. They do not establish the physical accuracy
of the fragmentation result or measure an isolated NLL correction.

## Numerical limits and remaining data work

Fixed-resolution domain changes are negligible at the displayed precision.
Doubled-grid responses reach about 0.3% of Born bins in photonic controls
and 0.6% in the coarser pointlike-HVP Z control. Some shifts are resolved;
there is no per-mille grid-convergence claim. Geometric-cutoff and angular
responses remain statistically noisy, with percent-level shifts and
comparable paired-iteration errors. Ten-iteration errors and cross-run
quadrature diagnostics are not certified numerical bounds. The PHOTOS
cutoff control measures PHOTOS, rather than its invariant analytic partner.

The auxiliary HVP input and 1.5–3 GeV matching-scale scans are small within
the declared pointlike prescription. They are not a complete assessment
of real hadronic-pair modelling. The original production-variance runs,
rejected grid/domain attempts and their explanations are preserved, while
the final report identifies exactly which runs enter each result.

Precision data use still needs a selected measurement, its energy/cuts and
prompt/extra-pair conventions, PDF/alpha_s and production-scale uncertainty,
and recoil treatment near the W Jacobian/low-qT regions. Finer grids and
better geometric sampling are needed for tighter numerical targets. Full
NLO EW, NNLO QCD and exact mixed corrections remain separate extensions.
Mixed strong/QED fragmentation kernels and fitted hadronic fragmentation
are not prerequisites to the present benchmark.

Fresh-run commands are in the [audit](../docs/qed_accuracy_audit.md).
`run_variance_ensemble.py` creates the ten-seed W− sample and its
`central_runs.json`; pass that file through `analyse_comparison.py --central-runs`.
For the actual retained campaign, also pass
`--domain-results results/qed_domain_v023` and
`--adapter-equivalence results/qed_adapter_regression_v023/report.json`.
Provenance and immutable raw data are retained alongside the figures.
