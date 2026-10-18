# NLO QCD production with finite W/Z decay matching

Version 0.16, 17 October 2026. The implemented
[contract](../docs/spin_decay_matching.md) is a leading-power unitary O(alpha)
decay, additively matched to fixed-alpha photonic LL. W± use the full charged-W
radiative current. Z/gamma* rates use its parity-even tensor and symmetric
pair cuts. H1 is not resummed in this increment; this is not a physical
multi-flavour full-NLL prediction.

## Inputs and observed rates

MCFM 10.3.1 NLO QCD, 13 TeV pp, NNPDF31_nnlo_as_0118 member 0 (the available
PDF validation input; this does not make the matrix elements NNLO), fixed
QCD and fragmentation scales MW=80.385 or MZ=91.1876 GeV. W is on shell;
Z/gamma* generation extends from 60 GeV to 13 TeV. Observed bare muons have
pT>25 GeV and |eta|<2.4, with no neutrino-pT/mT cut; Z pairs also require
60<m_mumu<120 GeV. No recombination, isolation or extra-lepton veto.

The FF comparison branch uses fixed alpha=1/137.035999084, m_mu=0.1056583755
GeV, mu0=m_mu, photonic non-singlet LL, 4096 log-loss bins over span 40, and
continuum MSbar D1 restoration. No pair channels or running are included.
Four independent angular-remainder samples are used per host entry, with
soft and 1-|c| exclusions 1e-5; these are numerical controls on a finite
remainder. Real and all mapped QCD dipoles share the original-real seed.

All rates are pb; errors below are host integration errors only. PHOTOS 3.56
uses the same physical alpha, interference on, W ME on and Z ME off, with a
1e-5 photon cutoff. Its exponentiated branch runs at physical alpha directly.
The first-order branch uses a statistically checked signed alpha/4 proposal,
as documented in the contract.

| Process | NLO QCD Born decay | Analytic additive NLO decay + LL | Paired exponentiated PHOTOS minus analytic |
| --- | ---: | ---: | ---: |
| W+ | 5488.56 ± 18.32 | 5371.51 ± 17.86 | -0.60 ± 5.40 |
| W- | 4187.43 ± 11.72 | 4104.02 ± 11.53 | +4.64 ± 4.45 |
| Z/gamma* | 750.15 ± 1.90 | 725.47 ± 1.79 | -1.71 ± 0.80 |

The directly accumulated **first-order PHOTOS minus analytic first-order**
differences are -4.88±7.15, +6.40±5.00 and -1.01±1.16 pb for W+, W- and Z,
respectively. These are separate estimators from the exponentiated comparisons.
The direct finite-hard additions are 52.34±2.23, 32.11±1.92 and 6.95±0.33 pb.

MCFM weights iterations independently for every histogram/bin. Consequently
its printed component means need not close by subtraction: use the directly
booked difference rows and their paired errors. There is no full cross-bin
covariance in this pilot. Neither agreement at this precision nor the roughly
2.1-sigma exponentiated Z difference establishes an NLL accuracy statement.
The branches retain different higher-order content and generator angular
approximations.

Nominal manifests and logs:

- `results/nlo_{wplus,wminus,z}_decay_order1/`
- `results/nlo_{wplus,wminus,z}_decay_exponentiated/`

Each records executable/core/adapter/PDF/configuration hashes and environment.
The complete source is frozen under `results/decay_source_v016/`; its core
source digest is `3acfea8f5028cf99d6f0345f6fe4c589cf814a6949573a7b77f556fcd527798a`.
The earlier failed first-order runs are retained under `*_spin_photonic/`
and `*_spin_photonic_v2/`, with their errors. They are excluded from results.

## Checks and numerical scope

All **50/50 CTest checks pass**, including the new polarised-current test and
massive W anchor (`results/decay_ctest.log`, 153 seconds). Checks include:

- Both W charges: photon Ward cancellation including the W graph, exact
  massive real energy spectra versus independent leading-power formulas,
  soft/collinear full-spin limits and all nine-component angular tomography.
- Z: independent existing vector-current trace, known energy spectrum,
  Hermitian rotation/contraction and mass-window geometric cancellation.
- Hard+FF second moments and fragmentation-scale cancellation at M/2,M,2M.
- Fixed-alpha photonic overlap removal and zero-coupling identity.
- Existing coupled kernels, evolution, endpoint/bin and reference checks.

The unpolarised W rest-frame anchor, with the same lepton pT/eta cuts, gives
C1=-14.05044061 from LP hard+FF and -14.05076175 from the exact massive real
current with unitary virtual normalisation. Their coefficient difference is
3.21e-4; the corresponding rates are 0.7666936852 and 0.7666904294. Tightening
both geometric exclusions by ten changes the geometry coefficient by
9.76e-6. These are tolerance/power/refinement observations, not certified
numerical bounds. See `results/spin_w_reference_checked.json`.

The optional PHOTOS proposal closure uses 50,000 paired Born points per
branch/process. Physical-alpha and signed alpha/4 first-order estimates differ
by -1.38 sigma for each W charge and +0.14 sigma for Z. This tests the sampling
estimator, not external matrix-element accuracy. The report and complete
per-point shifts are in `results/photos_linear_proposal/`.

Host controls use the `*_decay_zero`, `*_decay_slivers`,
`*_decay_grid8192` and `*_decay_alpha01` directories. The combined audit is
`results/decay_host_reference.json`, produced by `analyse_decay.py`. Its
separate-run comparisons use combined errors without cross-run covariance;
exact same-stream invariants are checked separately. Refinement shifts must
not be substituted for a complete numerical error budget.

## Reproduction

The host and PHOTOS builds are optional external dependencies. They are read
without modifying their installations. With the compatible local GCC core:

```sh
cmake -S . -B build-host -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 -DBUILD_TESTING=OFF
python3 adapters/mcfm/build_adapter.py /path/to/MCFM-10.3.1/build-gcc \
  --output build-mcfm-decay --photos-build /path/to/photos-build

python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --output results/my_wplus_order1 --adapter build-mcfm-decay/mcfm_leptonff \
  --threads 1 --photonic-ll --hard-samples 4 \
  --photos-order 1 --photos-linear-scale 4
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --output results/my_wplus_exponentiated --adapter build-mcfm-decay/mcfm_leptonff \
  --threads 1 --photonic-ll --hard-samples 4
```

Repeat for `wminus` and `z`. Every output directory must be new. Controls add
`--alpha 0`, `--hard-soft 1e-6 --hard-collinear 1e-6`, `--bins 8192`, or
`--dipole-alpha .1`. Keep the generator mode fixed when comparing a control.
`python3 adapters/mcfm/analyse_decay.py` audits the named production artifacts.
The core anchor is dependency-free: `./build/leptonff_spin`.

## Remaining phenomenology work

Dress H1 on all required hard parents for NLL and restore the corresponding
boundary-generated subleading towers; extend physical flavours/thresholds and
pair/charge counting; add physical hard-radiative pT, mT and pair-mass bins;
extend the independent polarised-amplitude grid; improve common-random sampling
and retain iteration/bin covariance; and perform QCD/PDF/fragmentation-scale
studies with a complete statistical/numerical budget. Only then use the result
as a controlled phenomenology prediction or for an experimental-data fit.


## Independent real amplitudes

`adapters/recola/check_spin.py` links a local RECOLA 1.5.1 and COLLIER 1.2.9
without changing either installation. Across 486 massless rotated/boosted
points and all three helicities, maximum relative squared-amplitude residuals
are 3.16e-14 (W-), 3.13e-14 (W+) and 5.59e-15 (charge-symmetrised Z). The
Z comparison isolates exactly the parity-even tensor required by these cuts.
See `results/recola_spin_checked/report.json` for the full provenance.

```sh
python3 adapters/recola/check_spin.py /path/to/recola1-1.5.1 \
  /path/to/COLLIER-1.2.9 --output results/my_recola_spin
```

## Component convergence and common-iteration audit

The final automated audit records **24 passed controls and one retained
component diagnostic**. Zero coupling, same-stream generator independence,
grid invariants and finite-remainder controls pass for W±/Z. Doubling the FF
grid changes the matched rate by +0.0547, +0.0815 and -0.04354 pb, respectively;
hard and PHOTOS columns are bitwise unchanged. These are measured shifts,
not residual error bounds.

The initial W- projected-hard row changes by 3.16 combined sigma under the
QCD dipole-alpha variation; the full matched row changes by only 0.72 sigma.
An independently seeded repetition with four times the initial calls gives
1.90 sigma for the projected row and 1.57 sigma for the full matched row.
The original flag is retained in `results/decay_host_reference.json` rather
than suppressed by changing its threshold.

An isolated MCFMStorage instrumentation then records the raw per-iteration
rate means and variances before the host applies different inverse-variance
weights to each histogram. Repeating the two high-statistics runs with this
instrumentation reproduces **all original printed histograms bitwise**.
Independent reconstruction from the raw records agrees within 4.5e-5 pb,
the printed precision. Common equal-iteration estimates restore component
algebra to 1.2e-10 pb. Their subtraction-parameter pulls are:

| W- component | Common-iteration pull |
| --- | ---: |
| Projected hard | -0.66 |
| Geometric hard | +0.41 |
| Complete hard addition | +0.24 |
| Matched rate | -0.18 |
| Paired PHOTOS minus matched | -0.66 |

At dipole alpha=1 the common-iteration matched estimate is **4121.97±13.71
pb**, and its paired PHOTOS difference is **-13.57±12.40 pb**. The corresponding
native independently weighted histogram gives 4109.76±6.62 pb and -3.90±2.60
pb on exactly those events. This difference exposes finite-sample weighting
and rare-event sensitivity. The smaller native errors in the initial rate
table must therefore not be promoted to a certified precision budget.

The raw audit supplies a full 19-rate covariance estimate from the ten
iteration chunks per QCD contribution, as well as the within-iteration
variances. It does not supply a controlled heavy-tail bound or cross-bin
covariance. See `results/wminus_iteration_comparison.json` and
`results/nlo_wminus_decay_audit{1,01}/iteration_rates.csv`. More effective
angular integration/sampling and stable joint estimators are still needed
for a precision phenomenology comparison.

Reproduce the instrumentation with `build_adapter.py --iteration-audit`,
then use `--seed 2089 --calls-factor 4` and dipole alpha 1 or .1 in the runner.
Audit the two new directories with:

```sh
python3 adapters/mcfm/analyse_iterations.py \
  results/nlo_wminus_decay_audit1 results/nlo_wminus_decay_audit01 \
  --output results/wminus_iteration_comparison.json
```

The additional adapter/reference source is frozen under
`results/decay_source_v016_audit/`; the core source digest is unchanged.
