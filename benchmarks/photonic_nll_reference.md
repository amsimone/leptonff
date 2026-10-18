# Photonic NLL implementation validation, v0.24

18 October 2026. This increment implements the explicitly
[unitary-conditioned fixed-alpha photonic NLL](../docs/photonic_nll.md)
mode. It does not replace the frozen v0.23 production comparisons.

| Validation | Result | Evidence |
| --- | --- | --- |
| Full core CTest regression, including new NLL matching/reference | 63/63 passed | `results/photonic_nll_ctest_v024.log` |
| Subsequent launcher and analysis CTests | 3/3 passed; one new test, two repeated analysis tests | `results/photonic_nll_pipeline_tests_v024.log` |
| Independent kernel/evolution reference | 55/55 checks passed | `results/photonic_nll_reference_v024.json` |
| Actual W+/W-/Z callback, LL/NLL routing and frozen LL regression | 39/39 checks passed | `results/photonic_nll_adapter_checked_v024.json` |
| Minimum-call NLO host, physical bins and grouped covariance | W+, W-, Z completed; 3/3 analysis checks passed | `results/photonic_nll_smoke_v024/report.json` |

The two CTest runs cover all 64 distinct configured tests. Independent
reference checks are also part of CTest, not additional independent evidence
to be counted twice. Four callback phase-space points per process reproduce
all 33 LL rates and 136 LL spectrum components bitwise against the frozen
v0.23 adapter/core. This is a callback regression, not a claim that we reran
and compared a complete LL production campaign.

The published source number moment is -0.7193116057745095. At physical
alpha and the nominal W fragmentation scale, the removed per-leg raw factor
is 0.9999871258728089; for Z it is 0.9999868811901362. These are the explicit
conditioning factors, not estimates of a physical inclusive weak correction.
At 4096 cells and span 40, the conditioned kernel and evolved number
residuals are about 1e-14, while the per-leg wrap/tail bounds on the added
NLL term are below 2e-17. These bounds do not include cell-discretisation,
hard-geometry cutoff or host Monte Carlo errors.

The independent continuum check refines from 1024 to 4096 to 65536 cells.
The coarse-grid differences in the second coefficient are retained; they
are not erased by endpoint fitting. Laboratory tests check the continuum
order-alpha result, the actual lattice second coefficient and cubic
remainders for W and Z. NLO corrections on H1/d1 and multiple subleading
insertions are excluded consistently with NLL counting.

The host smoke campaign uses `--call-scale .01`, four geometry samples and
4096 cells; **one W+ integration convergence warning is retained**. Its
cross sections and ratios must not be used for phenomenology. The report
writer and six-figure plot path were exercised in temporary output only.
No high-statistics photonic NLL plots were produced. The existing v0.23
plots are unchanged and still show photonic LL.

Production remains pending the reusable-record/binning work, finer numerical
controls and a statistical budget. No change to PHOTOS ME switches was made.
