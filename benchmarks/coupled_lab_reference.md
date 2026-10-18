# Coupled spin W/Z hard benchmark

Version 0.18 connects physical W± and parity-even, charge-symmetric Z/gamma*
decay currents to **one-family MSbar NLL** evolution on NLO QCD production.
Photon conversion, opposite-charge primary FFs, P1, running and d1 are all
included, with the common hard coupling at original M and continuum NLO
overlap restoration. See the [contract](../docs/coupled_lab_matching.md).
Physical flavours/thresholds, an on-shell input conversion, hard differential
spectra and a controlled phenomenology uncertainty budget remain open.

## Inputs and nominal rates

MCFM 10.3.1, 13 TeV pp, NNPDF31_nnlo_as_0118 member 0 as the available
validation PDF; matrix elements remain NLO QCD. Fixed QCD scales MW=80.385
and MZ=91.1876 GeV; default mu_frag equals that process's mass. The one-family
MSbar input is explicitly alpha(mu0)=1/137.035999084 at mu0=m_mu=0.1056583755
GeV, with two-loop one-family running. Its numerical equality to alpha(0)
is not a finite physical input conversion.

W is on shell. Off-shell Z/gamma* generation extends above the observed
120 GeV boundary. Muons are bare, pT>25 GeV, |eta|<2.4; Z observed pairs
satisfy 60<m_mumu<120 GeV. W has no neutrino-pT or mT cut. The coupled grid
has 4096 bins over log-loss span 40 and minimum observed fraction 1e-4.
Four geometric samples per QCD entry use common original-real seeds for all
mapped counterevents. Soft and 1-|c| geometric exclusions are 1e-5.

Nominal estimates use common equal-iteration means within fixed-call QCD
contributions, retaining zero iterations. Errors combine the raw
within-iteration variances. The ten-chunk covariance is retained as a
diagnostic; these are not certified rare-event, physics or cross-bin errors.

| Process | Coupled matched rate [pb] | Paired matched minus continuum order1 [pb] | Explicit a0² C2 [pb] |
| --- | ---: | ---: | ---: |
| W+ | 5382.935 ± 22.882 | -0.4064 ± 1.2251 | -0.5588 ± 1.4014 |
| W- | 4126.473 ± 36.486 | +0.6329 ± 0.6115 | +0.7166 ± 0.6289 |
| Z/gamma* | 725.2497 ± 1.8134 | +0.04448 ± 0.08757 | +0.04930 ± 0.09769 |

Component histograms have different native inverse-variance averages, so
native means/errors must not replace the common estimator or be subtracted
to make a paired uncertainty. Higher-order shifts are small compared with
the present integration errors; this sample does not establish their size
precisely. The large W- common error is retained.

The paired PHOTOS-minus-coupled differences are -2.415±6.315 pb (W+),
-38.047±32.235 pb (W-) and -2.1622±0.9828 pb (Z). PHOTOS 3.56 uses
exponentiation, fixed alpha(0), cutoff 1e-5, interference on, W ME on and
Z ME off. It has no pair channels. These are diagnostic differences between
distinct perturbative contents, not an aligned NLL scheme comparison.
The previous [photonic benchmark](photonic_dressing_reference.md) remains
the appropriate same-fixed-alpha comparison branch.

## Validation and controls

The final core suite passes 52/52 tests. The new checks cover:

- Complete photon projection plus geometry against the independently derived
  vector photon-pair coefficient for three pair-mass moments.
- Photon MSbar subtraction against the massive W/vector leading-power real
  spectra, and the unitary muon-plus-photon hard energy sum rule.
- Pointwise photon geometry against the independent Dalitz density.
- Literal signed hard-parent fraction sums, including conversions and wrong
  charges, on laboratory W/Z measurements.
- Common hard running against the independent RK4 moment evolution.
- Continuum O(alpha) overlap, zero coupling and explicit C2, with cubic
  higher-order remainders for both W charges and Z.
- W fragmentation-scale variations free of retained alpha² log(M/m_mu)
  terms, with a negative control that fails when photon parents are omitted.

The independent RECOLA 1.5.1 check retains all 486 helicity/rotation/boost
points: maximum relative differences are 3.16e-14 (W-), 3.13e-14 (W+) and
5.59e-15 (parity-even Z), in `results/recola_spin_coupled/report.json`.

All three nominal runs reproduce the previous continuum NLO, bare spin hard
and PHOTOS columns bitwise. W coupled FF columns reproduce the earlier
same-grid FF-only pilots; a new same-executable 4096-bin Z FF-only control
is used because the historical Z PHOTOS pilot used 2048 bins.

Runs are `results/nlo_{wplus,wminus,z}_coupled_{case}/`, with cases nominal,
zero, grid8192, alpha01, slivers, scale_half and scale_double. The last two
vary only mu_frag, keeping QCD scales and a_h(M) fixed. Required controls
check zero identities, invariant baseline columns, QCD dipole-alpha and
geometric refinements, and continuum NLO scale cancellation. Finite
matched scale/grid shifts are measured changes, not certified residual
error bounds. The combined report is `results/coupled_lab_host.json`;
nominal raw audits are in `results/coupled_nominal_iterations.json`.

## Reproduction and provenance

```sh
python3 adapters/mcfm/build_adapter.py /path/to/MCFM-10.3.1/build-gcc \
  --core-build build-host --output build-mcfm-coupled \
  --photos-build /path/to/photos-build --iteration-audit
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-coupled/mcfm_leptonff \
  --output results/my_wplus_coupled --threads 1 \
  --coupled-hard --hard-samples 4
python3 adapters/mcfm/run_coupled_controls.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-coupled/mcfm_leptonff
python3 adapters/mcfm/analyse_coupled.py
```

Repeat controls for W- and Z; every run directory must be new. The optional
iteration audit retains all 33 rate columns, including inactive-mode zeros.

The executable build source is frozen in `results/coupled_lab_source_v018/`,
core digest `a98e5026405774513b1637cd3aeb1b04c6d07c158afc60833b25e6f0d49de754`.
The final added W scale-cancellation/negative-control tests are frozen in
`results/coupled_lab_source_v018_final/`, core digest
`8b0bdbbe976fae115acbbb33f72282a9f5f9ac29927d03abd99cb4aa8bfd905a`.
These final changes affect tests only; the physics library hash is unchanged.
All 33 coupled bridge columns also pass finite/common-seed repeatability
checks over 1000 events for each W charge and Z, with the GCC host library;
see `results/coupled_smoke/report.json`. The new core test also passes when
compiled with the host GCC compiler.

Some early run manifests record the launch working-tree digest in their
source field; the embedded adapter build always retains the actual executable
source digest. The combined audit distinguishes those fields. Later runner
versions record executable and working-tree source hashes separately.
No older result is relabelled with the new test digest.

The zero-alpha control compares each branch to its own Born mapping:
massless for the analytic operator, massive for PHOTOS. The small retained
PHOTOS-minus-analytic Born difference is not assigned to radiation.

A further 16384-bin point is included for each process. Paired iteration
refinement diagnostics retain common QCD streams and verify every Born
iteration control exactly. Their ten-chunk errors are statistical diagnostics
only, distinct from the separate-run errors used in required checks.
The total-rate shifts are small, but the integration/grid resolution is not
yet sufficient to isolate the small higher-order contributions precisely.

## Completed numerical audit

All **30 host checks pass**, with no diagnostic flags. The maximum required
QCD-subtraction/geometric-refinement comparison is 1.42 combined standard
deviations. Continuum NLO fragmentation-scale cancellation passes for W±/Z.
The raw 33-column iteration identities close within 9.7e-11 pb.
There are 24 coupled host runs (eight cases per process) plus the same-grid
Z FF-only control. All runs completed without integration warnings.

The final core suite passes **52/52** in 138.43 seconds, in
`results/coupled_final_ctest.log`. The independent RECOLA and active-coupled
bridge smoke checks also pass.

The following are matched-rate changes from the 4096-bin nominal result.
Errors are the limited **paired iteration-chunk diagnostics**, not event-level
covariance errors or certified numerical bounds.

| Process | 8192 bins [pb] | 16384 bins [pb] | mu_frag=M/2 [pb] | mu_frag=2M [pb] |
| --- | ---: | ---: | ---: | ---: |
| W+ | +0.2906 ± 0.7605 | +1.1164 ± 1.0936 | +0.0346 ± 0.0717 | -0.0751 ± 0.0709 |
| W- | +0.1552 ± 0.2091 | +1.1633 ± 0.8833 | -0.0270 ± 0.0749 | +0.0003 ± 0.0794 |
| Z/gamma* | -0.0716 ± 0.0464 | -0.1181 ± 0.0548 | -0.0132 ± 0.0053 | +0.0055 ± 0.0053 |

The finest Z point is 0.118 pb below the coarse result, while the nominal
beyond-first-order shift is only +0.0445 pb. The W finest-grid changes are
about +1.1 pb with roughly 0.9–1.1 pb chunk scatter. No precise continuum
grid limit or sign of the small higher-order correction is established.
The much smaller fragmentation-scale band must not replace this unresolved
numerical/statistical assessment. Reducing the hard-estimator variance and
controlling the measurement/grid limit remain phenomenology tasks.

The final analysis/runner scripts and their hashes are archived in
`results/coupled_lab_analysis_v018/`. Per-run manifests retain the actual
build information and the source/working-tree distinction described above.
