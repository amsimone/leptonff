# Independent fixed-order photonic shape benchmark

The implemented normalized vector-current mass radiator has passed a new
external fixed-order test of its mass-logarithm content. McMule's strict
O(alpha^2) photonic FSR shape differs from our NLL expansion by a remainder
consistent with a constant in L=log(Q^2/m_mu^2), for each tested moment.
No logarithmic mismatch is resolved at this calculation's precision. This is
stronger physical evidence than checking our own formula against a second
implementation of the same formula.

The remainder is nonzero, changes with the moment, and is beyond NLL. It is
not added to the production kernel or interpreted as an absolute two-loop
width. The full resummed KKMC difference is still not completely explained.
The [full ingredient/evidence report](computation_report.md) places this test
alongside SANC, PHOTOS, PYTHIA, KKMC and the production interfaces.

## Process, source and perturbative definition

The reference is the pure vector-current process e+e- -> mu+mu- plus photons,
with radiation only from the muon line. The hard Q values are
10,30,60,91.1876,300,1000 GeV, and the muon mass is 0.1056583755 GeV.
Alpha is fixed at 1/137.035999084. ISR, IFI, vacuum polarization, real pairs,
fermion-loop running and weak corrections are absent. There are no angular or
lepton cuts. The physical massive support is used; no extra low-mass floor is
imposed. The analytic FF is its leading-power counterpart, so mass-power
effects remain explicit, particularly at the lowest Q.

The optional builder obtains the official
[McMule v0.8.0 source](https://gitlab.com/mule-tools/mcmule), pinned to
`de32e762ca44ab82345d4a511f7b64ab61ba9bd4`, with recursive pinned submodules,
including OpenLoops and COLLIER. The FKS-l subtraction and finite-mass treatment
are described in [Banerjee et al.](https://arxiv.org/abs/2007.01654).
The selected native pieces are:

| Piece | Use |
|---|---|
| `ee2mm0` | Born normalization check |
| `ee2mmFMM` | Inclusive muon-line one-loop plus integrated soft term |
| `ee2mmRMM` | Muon-line real emission, inclusive or weighted shape |
| `ee2mmRFMMMM` | Photonic muon-line real–virtual, with upstream subtraction |
| `ee2mmRRMMMM` | Muon-line double-real photons, with upstream subtraction |

No vacuum-polarization or pair pieces are included. No leptonff kernel enters
the external amplitudes. The Born denominator is the independently known
massive vector-current rate, checked by native integration:

    sigma0(alpha=1) = 4*pi/(3*s) * beta_mu/beta_e
                     *(1+2*m_e^2/s)*(1+2*m_mu^2/s).

The selected source sets internal alpha=1. After division by this Born rate,
relative NLO pieces are multiplied by physical alpha and NNLO pieces by
physical alpha^2. Comparing a coefficient in alpha with one in alpha/(2*pi)
requires the corresponding factor of 2*pi; the analysis does so explicitly.

## Why a double-virtual amplitude is unnecessary here

Let F(z)=z^(N-1), or theta(z>=cut), so F(Born)=1. Write

    sigma = sigma0+alpha*sigma1+alpha^2*sigma2
    M     = sigma0+alpha*M1+alpha^2*M2
    I1 = (M1-sigma1)/sigma0
    I2 = (M2-sigma2)/sigma0.

Then

    R = 1+alpha*I1+alpha^2*(I2-(sigma1/sigma0)*I1).

All no-real-photon terms, including double virtuals, cancel in the F-1 action
defining I2. Their absence is specific to this **normalized shape difference**.
It does not give an absolute NNLO rate. The inclusive NLO normalization
counterterm is indispensable and is integrated separately with disjoint seeds.

For stability the measurement factors F-1=(1-z)*[(F-1)/(1-z)]. The positive
loss 1-z guides VEGAS adaptation. For integer moments, the second factor is a
finite polynomial; it does not subtract nearly equal floating-point numbers.
The factorization is applied to **each native real and FKS counterevent**.
Soft subtraction is never replaced by an independently integrated singular
real-emission expression. The N=1 action vanishes exactly.

## Importance maps and build provenance

The native general phase-space sampler is retained as `--angular-map native`.
It is inefficient for high-Q photon–muon collinear configurations. The `fsr`
map rotates the same native phase space and resolves the mass-regulated
angular peaks logarithmically, with the explicit coordinate Jacobian. The
`fsr-direct` double-real option uses McMule's aligned two-photon map, samples
**both kinematic roots** and restores the actual incoming electron momenta.
Its photon energies, soft configurations and subtraction remain native.
These are integration coordinates, not changes to the matrix elements.

One million points per map at Q=3 GeV compare native and mapped three-/four-body
uncut Lorentz-invariant phase-space volumes, restoring the FKS photon energy
factors. Five kinematics checks and four volume checks pass. The largest volume
pull is 2.52; the two independent four-body importance maps agree. The physical
M_Z double-real results also agree between the mapped samplers. The weak native
pilot is retained as a control, not used to claim native-sampler ppm precision.

With gfortran 16 on this arm64 machine, the upstream generic two-body procedure
dispatch crashes while the direct call works. `prepare_source.py` inserts an
explicit dispatch to the original `psx2` routine for that case. The original
and prepared hashes are recorded. It changes no phase-space formula,
amplitude or subtraction. The virtual two-body piece does not enable a
real-emission near-to-soft branch; early failed initialization pilots are
excluded from all results.

The out-of-tree build sets `GIT_DIR`/`GIT_WORK_TREE` explicitly so upstream
version generation records McMule's commit rather than the enclosing leptonff
repository. The manifest includes recursive source hashes, driver/map hashes,
compiler information, the link command, dynamic-library hashes and the
Fortran-runtime hash. The runner verifies the recorded artifacts before and
after each job. It requires a complete finite native JSON result and no bad
phase-space points, because a Fortran STOP can return status zero.

The client uses its build tree's shared libraries. It is not described as a
self-contained sealed executable. Building a new driver gives it a content-
named executable, preserving an already running client's identity. Source and
linked-backend licensing remain separate from the original MIT client code.

## Completed runs and results

The main M_Z sample has two independent 20-million-requested-call real-shape
runs, four 24-million-call double-real runs, two 400,000-call real–virtual runs,
and separate Born/inclusive NLO integrations. Counts are requested VEGAS calls;
native stratification and soft counterevents change actual evaluations.
Five additional Q values use independent streams. Subtraction controls use
xi=0.05 and 0.5 against the central 0.2. A further control tightens softcut
1e-10 -> 1e-11, collcut 1e-14 -> 1e-15 and the real–virtual near-to-soft switch
1e-3 -> 1e-4 together. This last run is a combined technical control, not a
separate determination of each individual cutoff effect.

At Q=M_Z, define delta2 as the **strict alpha^2 contribution** from McMule
minus the alpha^2 expansion of our complete NLL-matched radiator. It is not a
difference between a fixed-order result and an unexpanded resummed curve.

| Moment N | O(alpha) difference [ppm] | delta2 [ppm] |
|---:|---:|---:|
| 2 | -0.650 +/- 0.642 | -28.073 +/- 0.719 |
| 3 | -1.438 +/- 1.116 | -17.041 +/- 1.012 |
| 5 | -2.938 +/- 1.983 | +5.566 +/- 1.513 |
| 10 | -5.470 +/- 3.792 | +47.224 +/- 2.634 |

Errors are independent-piece MC standard errors, including the measured
inclusive NLO counterterm's uncertainty. At M_Z the inclusive coefficient in
alpha is 0.238550+/-0.000397, consistent with its leading-power value
3/(4*pi)=0.238732. A normalization-only mismatch cannot describe the moment-
dependent signs in the last column.

The six-point constant fits to delta2 give:

| N | Constant [ppm] | chi-square / degrees of freedom |
|---:|---:|---:|
| 2 | -28.660 +/- 0.317 | 3.17 / 5 |
| 3 | -17.526 +/- 0.432 | 1.68 / 5 |
| 5 | +5.558 +/- 0.604 | 1.19 / 5 |
| 10 | +47.551 +/- 1.020 | 1.74 / 5 |

The remainder is compatible with a finite alpha^2 shape term. Fits allowing L
and L^2 are retained too; none resolves an additional logarithmic coefficient.
The N=2 all-point linear slope is 0.079+/-0.173 ppm per L. Excluding Q=10 GeV,
where mass powers matter most, it is 0.650+/-0.439 ppm per L. A linear-plus-
inverse-Q-squared sensitivity fit gives 0.771+/-0.517 ppm per L. A quadratic
fit has wider correlated coefficient errors. These numbers state the test's
precision and model dependence; they do not assert that every possible
logarithmic remainder is exactly zero. The largest-Q bins have larger errors.

For example the N=2 delta2 values are -28.624+/-0.410 ppm at Q=10 GeV,
-28.073+/-0.719 ppm at M_Z and -36.447+/-15.117 ppm at 1 TeV. The low-Q NLO
difference, +3.720+/-0.908 ppm, is a reminder that the massive and leading-power
calculations need not agree exactly. It is not hidden by resetting a mass or
fitting our first-order coefficient.

All three technical/subtraction controls agree. The largest pull across
the reported moments and cumulatives is 1.56, without treating those correlated
observables as independent global tests. Independent-seed scatter diagnostics
are retained per piece. MC errors come from native iteration scatter and are
not rigorous heavy-tail bounds.

![Fixed-order remainder versus the mass logarithm](../results/mcmule-fixed-order/comparison/fixed-order-log-scan.png)

The raw result also contains nine cumulative cuts. Precision near z=1 is much
worse than for the moments: at M_Z the z>=0.999 alpha^2 difference has an MC
error around 353 ppm. The campaign therefore validates the displayed moment
coefficients; it does **not** deliver a ppm-precision differential NNLO radiator
or certify every sharp endpoint bin.

## Interpretation for the earlier residual

This test materially reduces concern about the implemented photonic P1
mass-logarithm contribution for the inclusive vector observable. At the achieved
precision, it finds a finite second-order remainder rather than a demonstrated
NLL logarithmic failure. There is no evidence here for retuning the published
kernel or adding pair interference to make a photon-only comparison agree.

The finite remainder supplies useful physical information beyond our chosen
single-insertion NLL completion. Its signs are not reproduced by the simple
positive quadratic finite-factor variation for every moment. That variation
was a diagnostic, not a complete finite two-loop prediction or a bound.

The full resummed KKMC residual is larger and has mass dependence. One cannot
subtract this fixed-order number from that resummed residual and declare the
rest identified: higher-order dressings, exponentiation/residual conventions,
finite-mass effects and the CEEX support floor are still distinct. A full
finite NNLO matching function or controlled comparison of fixed-order KKMC
coefficients would be a subsequent calculation. No fitted constants have
been inserted into the production library.

This is also not an external NLL certification of all W angular/fiducial
observables, coupled pair-inclusive evolution or hadronic fragmentation. The
new evidence applies to the no-cut normalized vector-current mass observables
and the stated photonic content.

## Reproduce

Use a fresh result directory; the runners refuse to overwrite completed jobs.
Install NumPy/SciPy/matplotlib for analysis and meson/ninja/gfortran for the
optional external build. No PDF set, ROOT or production generator is needed.

```sh
python3 adapters/mcmule/build.py
python3 adapters/mcmule/check_maps.py --output results/my-mcmule/phase-space.json

python3 adapters/mcmule/run.py --output results/my-mcmule/normalization \
  --pieces ee2mm0 ee2mmFMM ee2mmRMM --mode inclusive \
  --calls 2000000 --iterations 10 --warm-calls 200000 --warm-iterations 8 \
  --seeds 2 --seed 921001 --jobs 4
python3 adapters/mcmule/run.py --output results/my-mcmule/real \
  --pieces ee2mmRMM --calls 2000000 --iterations 10 \
  --warm-calls 200000 --warm-iterations 8 --seeds 2 --seed 931001 --jobs 2
python3 adapters/mcmule/run.py --output results/my-mcmule/rr \
  --pieces ee2mmRRMMMM --calls 2000000 --iterations 12 \
  --warm-calls 200000 --warm-iterations 10 --seeds 4 --seed 941001 --jobs 4
python3 adapters/mcmule/run.py --output results/my-mcmule/rv \
  --pieces ee2mmRFMMMM --calls 50000 --iterations 8 \
  --warm-calls 10000 --warm-iterations 6 --seeds 2 --seed 951001 --jobs 2
python3 adapters/mcmule/analyse.py results/my-mcmule/normalization \
  results/my-mcmule/real results/my-mcmule/rr results/my-mcmule/rv \
  --output results/my-mcmule/central-report.json
python3 adapters/mcmule/campaign.py log-scan --output results/my-mcmule/log-scan
python3 adapters/mcmule/campaign.py controls --output results/my-mcmule/controls
python3 adapters/mcmule/summarize.py --root results/my-mcmule \
  --output results/my-mcmule/comparison
```

Build/run manifests record exact settings; different couplings, masses,
supports, importance maps, source builds, seeds or subtraction/near-to-soft
settings are guarded against inappropriate pooling. Source/driver changes
require a new build and distinct result groups. Log-scan analysis only permits
Q to change; technical controls keep physical settings fixed.

The local evidence is `results/mcmule-fixed-order/comparison/report.json`,
the phase-space check, PNG/PDF/SVG figures and native per-job JSON/log/checkpoint
records. The compact tracked snapshot is
[`benchmarks/mcmule_fixed_order.json`](../benchmarks/mcmule_fixed_order.json).
The complete existing research CTest suite passed 84/84 after repairing the
obsolete Python cache. The added optional-analysis CTest passes ten tests of
normalization, independent-product errors, pooling guards and logarithmic fits.
The Sphinx manual builds with warnings treated as errors. External integrations
remain optional, independently run evidence rather than CI dependencies.
