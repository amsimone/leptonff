# Observed W/Z spectra and covariance

Version 0.20 extends the NLO QCD host from fiducial rates to observed-momentum
histograms. The analytic fragmentation and PHOTOS measurements both act on
the final muon momenta. The old original-parent eta histograms are retained
as legacy controls; they are not used for the new radiative spectra.

## Measurement contract

The reference fiducial region is pT(mu)>25 GeV, |eta(mu)|<2.4, with no
neutrino or transverse-mass cut. Z pairs additionally require
60<m(mu mu)<120 GeV. The host retains the full original angular coverage
and, for Z, generates above the observed upper mass edge to allow migration.

| Process | Observable | Edges |
| --- | --- | --- |
| W± | observed muon pT [GeV] | 25, 30, 35, 40, 45, 60, infinity |
| W± | pT-based mT(mu,nu) [GeV] | 0, 40, 60, 70, 80, 90, 120, infinity |
| W± | observed absolute muon eta | 0, 0.6, 1.2, 1.8, 2.4 |
| Z/gamma* | observed pair mass [GeV] | 60, 70, 80, 85, 90, 95, 100, 110, 120 |
| Z/gamma* | minimum observed muon pT [GeV] | 25, 30, 35, 40, 50, infinity |
| Z/gamma* | maximum observed absolute muon eta | 0, 0.6, 1.2, 1.8, 2.4 |

Z entries count inclusive opposite-charge pairs from distinct original hard
parents. The pT and eta spectra above are symmetric pair observables, not
single-muon distributions or a leading-pair prescription. They respect the
charge symmetry required by the present parity-even neutral-current spin
tomography. Extra-pair vetoes and arbitrary experimental pair choices remain
outside this contract.

`LabSpectrum` expresses each bin as a difference of nested `LabSelection`
functionals. Threshold differences select (lower,upper]; eta differences
select [lower,upper). The fixed outer fiducial boundaries retain their strict
cuts, and disabled mT=0 retains the original no-cut definition. Endpoint
ownership is explicit even for numerical atoms. Every observable includes its
overflow, so its signed bin sum equals the fiducial rate without normalising
or clipping the result. JSON uses `null` for an infinite upper edge; overflow
integrals are reported without an invented finite density.

The mT convention is sqrt(2 pT_mu pT_nu (1-cos(delta_phi))). It also accepts
massive PHOTOS muons. The analytic fragmentation operators still separately
enforce massless parent momenta.

## Matching branches

The coupled branch uses the [one-family NLL contract](coupled_lab_matching.md)
and optionally [conservative fraction cells](conservative_fraction_cells.md).
Its eight columns per bin are Born, continuum order1, matched, retained
order2, matched-minus-order1, PHOTOS Born, PHOTOS post-FSR and the paired
PHOTOS-minus-matched difference. This PHOTOS comparison is diagnostic: the
running-coupling convention and pair channels differ.

The photonic branch uses [fixed-alpha LL with one finite insertion globally](photonic_lab_dressing.md).
Both tools then have fixed alpha and no pair production. This is an aligned
coupling/pair-content benchmark, not full NLL. PHOTOS exponentiation and the
analytic LL plus one finite insertion are different higher-order algorithms;
PHOTOS also retains muon masses, has the W matrix-element correction enabled
and the Z correction disabled. The auxiliary fourth analytic column is
explicitly labelled `ff2` and maps to the pre-existing FF-only LL expansion
coefficient. It is not a full dressed second-order coefficient. All other
columns have the same measurement meaning as in the coupled branch.

The 33 existing rate outputs retain their ABI and meanings. The additional
C entry point returns 8x17 physical bins. The original real momenta seed one
shared radiation stream for every mapped QCD contribution. All cumulative
measurements share each hard current, SO(3) rotation and collinear sector map.
PHOTOS generates one decay for the spectrum batch and measures every selection
on that decay. The legacy rate path is evaluated separately with the identical
seed, and its result must agree with the spectrum batch.

## Statistical contract

Each physical bin is booked as a separate one-bin MCFM histogram, with its
complete signed measurement weight. The host sums real and dipole entries
before forming squared weights. Native inverse-variance histogram means are
retained and independently reconstructed; subtraction of their independently
weighted means or errors is not used for physical comparisons.

The isolated adapter also records event cross-products for 144 channels:
eight total rates plus 136 spectrum columns. The event accumulator begins
and ends around one integrator call. Real contributions are accepted only at
the actual temporary-histogram commit site, after all dipoles. Failed real
points do not enter the products. Zero calls remain in the host call count.
All event sums and diagonal products are checked against the canonical host
histograms at every production iteration; disagreement stops the run.

For one iteration of N calls, with weights already including 1/N, the
covariance of its integrated estimates is

```text
C_ij = (N sum_events w_i w_j - I_i I_j) / (N-1).
```

Within each fixed-call QCD contribution, the analysis applies equal weights
to all production iterations, including zero iterations. It then adds the
independent contribution covariances. Every spectrum's summed covariance is
checked against every total-rate covariance entry. Normalised bin fractions
use the joint bin/total covariance, including the denominator fluctuation.
The paired PHOTOS difference has directly accumulated statistics.
Full normalised-shape covariance matrices are also exported and checked
against both the propagated diagonal errors and the sum-to-one constraint.
The combined matrix has exact linear dependencies; fits must select an
independent set of bins/components or account for those constraints.

Ten-iteration fluctuation covariances are retained separately as diagnostics.
Full event cross-products remove the previous cross-bin information gap;
they do not certify rare-tail convergence or make a ten-iteration sample a
precision error theorem. A control-variate estimator can give a noisy signed
estimate even in a kinematically forbidden tail whose expectation is zero.
No such estimates are clipped. The covariance implementation currently
requires one thread, central-PDF histograms and the audited NLO W/Z host.
Separate PDF-member runs are supported.

## Running and analysing

Build with the compatible GCC core and the isolated host instrumentation:

```sh
python3 adapters/mcfm/build_adapter.py /path/to/MCFM-10.3.1/build-gcc \
  --core-build build-host --output build-mcfm-spectra \
  --photos-build /path/to/photos-build --iteration-audit
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-spectra/mcfm_leptonff --threads 1 \
  --output results/my_wplus_spectra --spectra \
  --coupled-hard --fraction-cells --hard-samples 4
python3 adapters/mcfm/analyse_spectra.py results/my_wplus_spectra \
  --output results/my_wplus_spectra.json
python3 adapters/mcfm/plot_spectra.py results/my_wplus_spectra.json \
  --output results/my_wplus_spectra
```

For the aligned photonic branch, replace `--coupled-hard --fraction-cells`
with `--photonic-ll --dress-hard`. Repeat for `wminus` and `z`.
Plots are standalone PNG/SVG files; only plotting requires Matplotlib/NumPy.

`run_spectrum_controls.py` runs all three processes with explicit call budgets,
zero coupling, grid/cutoff refinement, an independent seed, three fragmentation
scales and seven QCD scale points. `--model photonic` selects the aligned
branch. The runner refuses to overwrite runs and checks completed-run
provenance when resuming. `--call-scale .05` is a reduced validation budget;
the default single-run budget is one. QCD renormalisation/factorisation scale
factors and individual PDF members can also be specified in `run_reference.py`.

The coupled campaign analyser reports the QCD envelope, fragmentation envelope,
grid and geometric-cutoff sensitivities, seed shifts and integration errors as
separate components. It does not combine unvalidated sensitivities into a
certified uncertainty. PDF replicas, physical flavour/input matching and
precision recoil treatment remain separate phenomenology gates; see the
[operational summary](operational_summary.md).
