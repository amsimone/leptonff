# NLO QCD W/Z host and PHOTOS pilot

Version 0.15.0 · 17 October 2026 · MCFM 10.3.1 / GCC 16.2 · 53-bit Real

The [laboratory contract](../docs/lab_host_matching.md) defines the new action.
This is an operational **NLO QCD plus Born-parent FF-only** benchmark with a
same-production PHOTOS comparison. The process-specific finite O(alpha) W/Z
hard-decay coefficients are absent. These numbers are implementation evidence,
not fully matched physical NLL predictions or measurements of an NLL effect.

## Inputs and comparison

All runs use 13 TeV pp and muon pT>25 GeV, |eta|<2.4, with no neutrino-pT or
mT cut. Z/gamma* additionally requires observed 60<m_mumu<120 GeV. W is in the
pole/zero-width approximation at 80.385 GeV; Z/gamma* uses the host's off-shell
process with MZ=91.1876 GeV. The production matrix elements retain spin
correlations and NLO QCD recoil. MCFM EW corrections are disabled.

The available `NNPDF31_nnlo_as_0118`, member 0, is used consistently for host
validation; the matrix elements are NLO QCD. This does not upgrade the production
accuracy. QCD renormalisation and PDF factorisation scales are fixed to MW/MZ.
FFs use the same central fragmentation scales, mu0=m_mu=0.1056583755 GeV,
alpha_MSbar(mu0)=1/137.035999084, the closed one-family coupled NLL evolution,
log span 40 and a supported minimum fraction of 1e-4.

PHOTOS 3.56 uses exponentiation and interference, W ME correction on and Z ME
correction off, alpha(0)=1/137.035999084 and photon cutoff 1e-5. It acts on the
same signed contributions before cuts, with a common photon random stream for
each real/dipole group. Its physical-muon-mass Born projection is recorded
separately. No QCD shower or duplicate FSR is applied.

## Pilot rates

All entries below are **pb**, with MCFM integration errors only. W rows are
tagged by the hard W/primary-lepton charge. Baseline FF grids have 4096 points
for W and 2048 for Z. The [combined artifact](../results/nlo_host_reference.json)
contains the detailed run manifests, inputs, configurations and hashes.
The [frozen implementation snapshot](../results/nlo_host_source_v015/CMakeLists.txt)
preserves the core and adapter sources before further hard-matching work.

| Process | Native-measurement Born parent | FF-only O(alpha) restored resummation | PHOTOS post-FSR |
| --- | --- | --- | --- |
| W+ → mu+ nu | 5488.6 ± 18.3 | 5314.3 ± 16.7 | 5366.8 ± 18.3 |
| W− → mu− nubar | 4187.4 ± 11.7 | 4070.0 ± 10.8 | 4113.9 ± 11.0 |
| Z/gamma* → mu mu | 750.15 ± 1.90 | 718.19 ± 1.71 | 723.54 ± 1.88 |

The directly accumulated PHOTOS-minus-FF estimates are 53.57±4.89 pb, 36.16±4.20
pb and 5.694±0.763 pb, respectively. These are comparison-pipeline diagnostics.
The FF branch lacks finite hard terms and has different scheme/pair content;
the differences cannot be interpreted as an NLL improvement or generator error.

MCFM weights iterations independently per histogram/bin. Thus the printed
central rates do not close algebraically with the directly accumulated
difference estimators. Use the latter for paired uncertainties. Full cross-bin
covariance and high-statistics bias/convergence studies remain open.

## Verification

- All **48 CTest checks pass** in the [Release test log](../results/lab_fragmentation_ctest.log).
  The new test compares the W action with literal signed lattice sums and
  independent continuum distribution quadrature, tests neutrino/mT kinematics,
  checks Z's full same/opposite tensor and C2, and tests zero-coupling and cubic
  perturbative remainders.
- With zero QED coupling, W and Z Born/expanded/resummed/restored FF rate
  columns agree exactly at printed precision. Native host agreement is 0.16,
  −0.89 and 1.54 combined standard errors for W+, W− and Z, respectively.
  This is a statistical interface check, not a per-mille accuracy certificate.
- Changing the QCD dipole alpha from 1 to 0.1 keeps the tested W+/W−/Z Born,
  FF-restored and PHOTOS rates and paired differences within 1.21 combined
  standard errors.
- Doubling the W+ grid 4096→8192 changes its FF-restored rate by +0.1933 pb;
  W− 4096→8192 changes by +0.3119 pb; Z 2048→4096 changes by +0.1745 pb.
  The native-production and PHOTOS rate
  columns remain identical under those grid changes, verifying common production
  and photon streams. These are measured shifts, not error bounds.
  Differential-bin refinement budgets remain open.
- Optional PHOTOS smoke tests check both W charges and Z for active radiation,
  momentum conservation and common-seed repeatability. The
  [build/test log](../results/photos_adapter_build_and_test.log) records them.
  MCFM-to-massless energy repairs are monitored and bounded; sampled maxima
  are recorded in each run log.

Seeded pseudorandom host integration is used for the accepted comparison runs.
Earlier Sobol pilots had large final W integration chi-squared and are not the
comparison central points. Runs marked `INVALID_HISTOGRAM_INDEX.txt`, failed
manifests, and native-only exploratory runs are not physics benchmarks.

## Reproduction

The optional adapter reuses a local compatible host build without changing it:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 4
ctest --test-dir build --output-on-failure --parallel 4
cmake -S . -B build-host -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-16 -DBUILD_TESTING=OFF
python3 adapters/mcfm/build_adapter.py \
  /Users/simone.amoroso/Work/generators/MCFM-10.3.1/build-gcc
python3 adapters/mcfm/build_adapter.py \
  /Users/simone.amoroso/Work/generators/MCFM-10.3.1/build-gcc \
  --output build-mcfm-photos \
  --photos-build /Users/simone.amoroso/Work/generators/w-weak-reweight/build/photos
python3 adapters/mcfm/run_reference.py \
  /Users/simone.amoroso/Work/generators/MCFM-10.3.1 wplus \
  --output results/new_wplus_photos --threads 1 --bins 4096 \
  --adapter build-mcfm-photos/mcfm_leptonff
```

Use `wminus` or `z` for the other processes, `--native` for native cuts,
`--alpha 0` for identity, `--dipole-alpha .1` for the subtraction check, and
double `--bins` for grid checks. Core-only runs use the default adapter and
need no PHOTOS dependency. The runner refuses to overwrite a result directory,
rejects bypassed callbacks/non-finite histograms and records final integration
chi-squared. `analyse_reference.py` collects the named reference runs and checks
the identity/control evidence; it does not assign a physical accuracy label.

Next: align photonic/flavour/coupling and decay-width conventions, implement the
spin-correlated finite W hard match and then the Z counterpart, extend observed
pT/mass/mT measurements, and quantify the physics and numerical uncertainty
budget. NNLO QCD follows the validated NLO W/Z programme.
