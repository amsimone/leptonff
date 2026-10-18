# Z mass radiator: validation and uncertainty diagnostics

This study validates the **stated normalized, fixed-alpha photonic prescription**
separately from its comparison with an independent physical photon generator.
It stays on Born QCD. No W production is launched by this study.

The implementation and retained logarithmic coefficients pass the independent
audit. KKMC comparisons resolve small residual differences; those are reported,
not converted into an automatic claim of physical NLL validation or a total
theoretical error.

The subsequent [residual investigation](z_radiator_residuals.md) adds coupling
and mass controls. It resolves a mass-dependent remainder, records the legacy
1 GeV CEEX support floor, and explains why fixed-order coefficient separation
is needed before assigning the discrepancy to a particular missing term.

## Complete observable, not FF-only coefficients

At fixed parent invariant mass Q, define z=m_mumu^2/Q^2 and the Mellin moment
R_N = integral dz z^(N-1) R(z). Let

    a = alpha/(2*pi)
    L = log(Q^2/m_mu^2)
    Lf = log(mu_frag^2/m_mu^2)
    ell = L-Lf
    P_N = integral dz p(z)*(z^(N-1)-1), p(z)=(1+z^2)/(1-z)
    C_N = integral dz p(z)*(log(z)-1)*(z^(N-1)-1)
    K_N = integral dz G_regular(z)*(z^(N-1)-1).

C_N is the sum of the one-loop initial boundary and per-leg projected hard
coefficient at mu_frag=Q. Their endpoint logarithms cancel in this sum. The
implemented single-insertion continuum expression is

    U_N = exp(2*a*Lf*P_N)
    R_N^LL = U_N * [1+2*a*(C_N+ell*P_N)]
    R_N^NLL = R_N^LL + 2*a^2*Lf*U_N*K_N.

The factor two includes both measured lepton legs. Its complete second-order
coefficient is

    r1_N = 2*(L*P_N+C_N)
    r2_N^LL = 2*Lf^2*P_N^2 + 4*Lf*P_N*(C_N+ell*P_N)
    r2_N^NLL = r2_N^LL + 2*Lf*K_N.

At mu_frag=Q the retained L^2 and L coefficients are respectively
2*P_N^2 and 4*P_N*C_N+2*K_N. This includes the hard coefficient; the library's
historical FF-only auxiliary second coefficient is not used as a full
observable reference.

For integer N, P_N and C_N follow from rational polynomial integrals:

    p(z)*(z^(N-1)-1) = -(1+z^2)*sum(z^j, j=0..N-2).

For example P_2=-4/3 and C_2=22/9, independently reproducing the normalized
first invariant-mass moment r1_2=44/9-8*L/3.

The independent Python checker uses these rational expressions, a separately
implemented Gauss integral and SciPy's dilogarithm for G. Direct double-plus
quadratures check cumulative distributions z>=z_min, independently of FFT
transport and the moment factorization. The native diagnostic evaluates the
same piecewise constant log-loss density as the production mass-bin operator,
including its endpoint atom.

The audit scans Q=75,91.1876,120 GeV, muon-mass factors 1/2,1,2,
fragmentation-scale factors 1/2,1,2, and alpha factors 1,1/8,1/16,1/32.
It checks moments N=1,2,3,5,10 and cumulatives at z_min=0.25,0.5,0.9,0.99,0.999.
Independent 65,536/131,072-cell calculations control numerical differences.
The complete second coefficient is extracted by a small-alpha extrapolation;
the remaining term is checked to be cubic. A mass-log fit checks the L^2,
L and constant pieces separately. Compensated moment summation is necessary
to avoid magnifying roundoff in this extraction.

## Normalization audit and its limit

The time-like photonic diagram class G is taken from
[Arbuzov–Melnikov, eq. (15)](https://arxiv.org/abs/hep-ph/0205172).
Its independently integrated first moment is

    G_1 = -13/4 + 3*zeta(2) - 2*zeta(3) = -0.719311605774...

The source is left intact; the normalized shape uses K=G-G_1*delta.
At two-leg level the raw source-sector inclusive factor is
1+2*a^2*Lf*G_1 to the retained order. Expanding division by that factor gives
the displayed K_N contribution. Exact unexpanded division would also generate
terms beyond NLL.

This verifies the normalization algebra of the chosen source sector. It does
not derive an absolute two-loop partial width or establish that every
generator's definition of photon-only radiation has the same sector
bookkeeping. In particular, imposing R_1=1 does not validate the shape.
The cancellation involving photonic and pair-interference classes in
[Arbuzov, section 5.2](https://arxiv.org/html/hep-ph/0206036v2) must not be
replaced by a claim that G is an independently conserving valence kernel.

## A controlled matching-prescription variation

The optional mass-only `quadratic_finite` prescription adds

    delta R_N = 2*a^2*U_N*(C_N+ell*P_N)^2.

It is the quadratic term of a finite-factor exponentiation, with the rest
of the default calculation retained. It changes a finite alpha^2 term,
preserves first order and the NLL mass-log coefficients, and preserves
normalization. It is one NNLL completion, not a full NNLL calculation or
a proof that all unknown terms are bounded by its size.

The same term is implemented as causal, zero-padded convolutions for mass bins.
The scale dependence of the quadratic term is retained consistently. No
P1-on-H1 or P1-on-boundary products are added.

## Independent KKMC comparison

The optional adapter builds the official Fortran KKMCee v4.32.01 source at
commit `527ed9d6f8c0080a256e77ae0d7adbf4c5207eea`. It needs C++ and gfortran,
but neither ROOT nor EW lookup tables. No upstream amplitude is modified.
Source and driver hashes and the build command are retained.

At fixed Q it generates muon pairs with fixed on-shell alpha and no ISR,
initial-final interference, weak corrections, hadronization, real pairs or
the optional virtual-pair correction. These choices are checked in the
upstream switches and order-weight assignments. The normalized observables
use each perturbative order's own inclusive weight.
The pinned defaults also impose MinMassCEEX=1 GeV: the CEEX denominator
excludes the small lower-mass tail. The leading-power first-order estimate
of that omitted fraction at Q=M_Z is 0.70 ppm, below the observed residual.
It remains a support convention to align for a definitive comparison; see
the residual investigation for the estimate's limits.

Correlated weights 72/73 are EEX residual orders 1/2; 251/252/253 are CEEX
orders 0/1/2 with interference off and without the pair-added weight 263.
CEEX2 contains second-order subleading QED residuals; see the
[KKMCee paper](https://arxiv.org/abs/2204.11949).
**EEX1/CEEX1 are exponentiated predictions, not fixed O(alpha), and these
order labels are not synonyms for leptonff LL/NLL.**

Two upstream restrictions are explicit:

* CEEX requires removal of photons below its technical soft threshold from
  the visible event record. Their recoil remains in the muons. The adapter
  bounds the missing soft four-vector instead of demanding exact closure of
  that incomplete record.
* With both ISR and FSR off, CEEX partitioning is unavailable in this release.
  The strict radiation-off control therefore exercises EEX; CEEX instead has
  a separate near-zero-coupling control.

Legacy Fortran STOP can return status zero. The runner requires a complete
result file too, so such a stop cannot masquerade as a successful reference.
Each campaign snapshots its executable and defaults before generating events.

## Results and interpretation

The sealed full audit passes **4,623 checks**. Targeted CTests cover the
mass response, photonic kernel/dressing, covariance analysis and the
independent complete-radiator checker. The small CI version uses two coarser
grids and the same physical tests when NumPy/SciPy are available.

KKMC has 200 million weighted trials at each of the three Q values, a
100-million-trial cutoff control at Q=M_Z, and 20 million trials at each of
alpha/8 and alpha/16 for the first-order limit. Radiation-off and near-zero
coupling controls have 10,000 trials each. All within-observable cross-order
covariances are retained. Each central and cutoff result combines ten
independent seeds. Moments and nested cumulatives are correlated;
their individual pulls are not a global chi-square.

The tenfold soft-cutoff change has maximum pull about 1.02 over the monitored
observables/orders. The coupling extrapolation agrees with the independently
known first-order moments at the available MC precision; it is not a
replacement for the much sharper massive SANC first-order check.

At Q=M_Z, the N=2 CEEX2 moment is 0.96447111 +/- 0.00000713, compared with
0.96455716 from leptonff NLL: a relative difference
-0.008922% +/- 0.000739% (MC only). The largest absolute difference among the
tested moments is 1.31e-4. Some cumulative residuals are also resolved:
z>=0.25 differs by -0.007801% +/- 0.000481%, while z>=0.999 differs by
+0.029503% +/- 0.004736%. The tenfold statistics increase reduces MC errors
by factors 3.16–3.41. The original sample versus the independent added sample
has maximum absolute pull 1.88. Seed-versus-rest pulls have maximum 3.22
across 390 correlated, overlapping checks (RMS 0.94).

A separate million/two-million-cell profile controls the library's endpoint
integration. Its refinement shift at z>=0.999 is 0.001996% relative to NLL,
about 0.42 times the current MC error; it is displayed separately and is an
empirical numerical diagnostic, not a rigorous error bound. The residuals
persist with the increased statistics; their physical interpretation remains
separate from their MC significance.

These differences support agreement at the few-1e-4 level for these
observables, but not equality at 1e-4 precision. They are not uniformly
bracketed by the current scale and matching diagnostics. Missing finite
second-order terms, different higher-order soft/collinear completions and
the sector/normalization translation need to be distinguished before
attributing a residual to one cause. The NLL increment is not automatically
a better approximation to every generator observable.

For the existing pp Born-QCD spectrum integrated over 90–92 GeV:

| Diagnostic | Relative change from default NLL |
| --- | ---: |
| mu_frag=M_Z/2 | +0.01169% |
| mu_frag=2*M_Z | -0.03375% |
| quadratic finite completion, mu_frag=M_Z | +0.01258% |

The additional NLL correction relative to LL is about +0.01381% there.
Thus a legitimate beyond-NLL completion can be comparable to that small
increment. This is a reason to retain the distinction between resolving an
increment numerically and establishing its physical accuracy.

No total theoretical uncertainty is assigned. Pair/running-alpha, ISR/IFI,
weak corrections and improved QCD production remain excluded physical
contributions, not entries implicitly covered by the QED scale envelope.
Finite-mass higher orders and an absolute branching-fraction normalization
also require their own treatment for phenomenology.

## Reproduction and artifacts

Install the optional analysis dependencies with
`python3 -m pip install -r adapters/validation/requirements.txt`.

```sh
cmake --build build-architecture --target leptonff_mass_audit
python3 adapters/validation/check_mass_radiator.py \
  --client build-architecture/apps/leptonff_mass_audit \
  --output results/my-radiator-audit

python3 adapters/kkmcee/build.py
python3 adapters/kkmcee/run.py --output results/my-kkmc-fixed-q \
  --Q 75 91.1876 120 --trials 20000000 --jobs 2
```

The KKMC runner also accepts `--cutoff 1e-6`, `--alpha-factor`,
`--mass` and `--fsr-off` for separate controls. The analytic Z runner accepts
`--mass-matching quadratic_finite` together with its fragmentation-scale scan.
Do not average these prescriptions or numerical refinements as independent
statistics.

### Independent seeds and increased statistics

Use `--seeds` to generate independent chunks in isolated directories, for example
`--trials 20000000 --seeds 2000001 3000001 4000001 --jobs 3`.
Each parent-mass index adds 100003 to its base seed; the runner rejects
duplicate resolved seeds. Retain the same executable and defaults snapshots
when extending an existing campaign.

The comparison accepts several campaign directories after `--central` or any
control option. It rejects overlapping seeds, different physics/cutoff settings,
or different executable/defaults hashes within a pooled sample. It combines raw
trial-weighted means and covariance matrices before forming normalized ratios.
For chunk sizes n_i, means m_i and covariances C_i **of the means**, it uses

    n = sum_i n_i
    m = sum_i n_i*m_i/n
    C = [sum_i n_i*(n_i-1)*C_i
         + sum_i n_i*(m_i-m)*(m_i-m)^T] / [n*(n-1)].

This is the covariance obtained by concatenating the independent trial samples;
it retains both within-chunk and between-chunk scatter and cross-order weight
correlations. Separate seed-versus-rest CEEX2 pulls check reproducibility.
Those tests overlap, and the observables are correlated, so their maximum pull
is a diagnostic rather than a global goodness-of-fit probability. Numerical
refinements and changes of physical prescription are never pooled as MC samples.

Local results are in `results/z-radiator-audit/`:

* `complete-sealed/`: independent audit inputs, executable snapshot and report.
* `comparison-high-stat/`: fixed-Q moment/cumulative figures, KKMC controls,
  seed checks and precision gain relative to the original 20-million-trial samples.
* `uncertainty-high-stat/`: folded Z scale/matching figures and updated diagnostics.
* `kkmc-central-high-stat/`, `kkmc-cutoff-high-stat/`: additional independent
  chunks, combined with the original `kkmc-central/` and `kkmc-cutoff/` samples.
* `comparison-final/`, `uncertainty-final/`: archived original-statistics reports.
* `fixed-q-profile/`: fine endpoint-grid controls.

The small tracked snapshot is
[`benchmarks/z_radiator_validation.json`](../benchmarks/z_radiator_validation.json).
Generated samples, plots and external sources remain untracked.

Remaining scientific validation: classify the resolved external residuals
through a physical second-order benchmark with identical sector conventions;
then fold a sufficiently controlled independent radiator with the same pp
production. The current fixed-Q moment differences must not be copied into
an mll-dependent phenomenological error band.
