# NLO+NLL QCD with photonic NLL decay: inclusive dimuon mass

The first combined example uses angular-inclusive Z/gamma* production from
DYTurbo at NLO+NLL and the leptonff fixed-alpha photonic decay prescription.
It compares QED off, the finite order-alpha decay correction, and matched
photonic NLL on one common production spectrum. There are no laboratory muon
pT or eta cuts. This is a restricted, useful first coupling; the general
native W differential interface is described [separately](dyturbo_w_decay.md);
general fiducial Z coupling remains open.

For a mass-only measurement the integrated decay response is rotationally
invariant. The trace of the finite angular remainder vanishes at each radiation
point; `inclusive_mass_response` tests this independently for both a bin
measurement and a smooth evolved measurement. The inclusive projected hard
coefficient and its evolved convolution therefore supply the finite matching
without needing the production spin density. QCD real and mapped subtraction
events share the boson invariant mass, so the same mass response applies to
the resummed, asymptotic and finite production contributions.

The implementation uses a bin-integrated production spectrum with a midpoint
mass-migration operator:

\[
 \sigma_i \simeq \sum_j \sigma_j^{\mathrm{NLO+NLL\ QCD}}
 R_i(M_j),\qquad M_j=(M_{j,\min}+M_{j,\max})/2.
\]

Each parent-mass bin migrates into all lower observed mass bins. Reusing the
production spectrum for all three decay prescriptions makes their difference
depend on the decay alone. Midpoint evaluation is a numerical approximation;
0.5 and 0.25 GeV production grids and doubled fragmentation grids quantify it.
This does not yet provide a continuous native-integrand decay callback.

Settings: pp at 13 TeV, the same NNPDF31_nnlo_as_0118 central member and native
DYTurbo EW/scale/subtraction/resummation settings as the production benchmark;
66 < parent mass < 116 GeV. Decay alpha is fixed at 1/137.035999084,
mu0 = m_mu = 0.1056583755 GeV, fragmentation scale 91.1876 GeV. The photonic
NLL kernel is unitary-conditioned and has one subleading insertion globally,
including finite order-alpha matching. Real pairs and fermion-loop running
are absent in this model. Observed mass bins are 1 GeV, from 50 to 116 GeV.
The complement below 50 GeV is assigned by the unitary normalization convention;
that sum is an arithmetic closure check, not an independent normalization test.

The positive pT support used internally by the generic decay table is a dummy
rest-frame bound on momentum fractions of 1e-4. It cannot restrict these visible
bins: m_pair >= 50 GeV and parent mass <= 116 GeV imply each fraction > 0.18.
It does not introduce a laboratory fiducial cut.

Build the optional example against an installed leptonff library:

```sh
cmake -S adapters/dyturbo -B build-dyturbo-interface \
  -DDYTURBO_ROOT=/path/to/dyturbo/prefix \
  -DCMAKE_PREFIX_PATH=/path/to/leptonff/prefix \
  -DLEPTONFF_DYTURBO_MASS_EXAMPLE=ON
cmake --build build-dyturbo-interface --parallel 4
python3 adapters/dyturbo/run.py --mode nlo-nll \
  --mass-bin-width .25 --accuracy .002 --output results/z-production-fine \
  --pdf-data /path/to/PDF/grids
# Repeat with .5 for the production-grid control.
python3 adapters/dyturbo/compare_mass.py \
  --production results/z-production-fine/nlo-nll/native.json \
  --coarse-production results/z-production-coarse/nlo-nll/native.json \
  --decay-executable build-dyturbo-interface/leptonff_dyturbo_mass \
  --output results/z-qed-effect --grid-bins 8192
```

The comparison script needs NumPy and Matplotlib. It writes spectra, numerical
controls, a JSON report and PNG/SVG plots. It retains signed corrections and
does not infer an event covariance from quadrature bins. Propagated numerical
error estimates are not complete QCD/QED theory uncertainties. Production
outside 66..116 GeV is excluded, so feed-in from larger parent masses is absent.
An experimental fiducial comparison needs a wider production domain and the
angular coupling, followed by PDF/scale/cutoff and numerical precision studies.

## First numerical comparison

The [benchmark record](../benchmarks/dyturbo_qed_mass.json) uses 200 production
mass bins and 65,536 fragmentation nodes. In the observed 90..92 GeV window:

| Decay prescription | Cross section [pb] |
| --- | ---: |
| QED off | 811.455 |
| Through order alpha | 697.443 |
| Matched photonic NLL | 703.723 |

The full photonic correction lowers this peak window by 13.28%; higher orders
raise it by 0.90% relative to the order-alpha result. Radiation enhances the
lower-mass tail. For the 66..116 GeV observed window the corresponding rates
are 1878.492, 1835.587 and 1834.657 pb. The unitary fully inclusive convention
redistributes probability into lower masses; a restricted window loses events.

The peak NLL/QED-off ratio changes by 0.055 percentage points between 0.5 and
0.25 GeV production grids, and 0.003 percentage points between 32,768 and
65,536 fragmentation nodes. Across visible bins the largest changes are about
0.41 and 0.10 percentage points respectively. These are refinement differences,
not rigorous uncertainty bounds. An explicit alpha=0 run closes bin by bin.
The local plot/data/report are under
`results/dyturbo-qed-effect/comparison-final/`.

## Fine-bin comparison with PHOTOS

`mass_precision.cpp` adds a conservative quartic reconstruction of the native
bin-integrated mass spectrum and integrates the decay within each source bin.
Integration intervals split at every observed mass edge. A zero-padded FFT
reduces the two-leg lattice action to a cumulative total-loss table, with one
finite insertion globally. Continuum first-order matching and the Born atom
are retained exactly. Log-loss cells are averaged for interpolation; their
spacing must be refined. This is a private adapter reduction, checked against
the public two-leg operator, at zero alpha and in its small-alpha expansion.

The source-grid error propagation includes all five bins in each conservative
reconstruction stencil, rather than treating subnodes as independent production
samples. Native quadrature errors are numerical estimates, not random event
variances or QCD theory errors.

The optional `leptonff_dyturbo_mass_precision` executable is built with
`LEPTONFF_DYTURBO_DECAY_CLIENT=ON` and the external PHOTOS inputs described in
the W guide. It takes a four-column source table (`lower upper integral error`,
in GeV/fb), a new output JSON path, fragmentation nodes, quadrature nodes,
PHOTOS trials, seed and alpha:

```sh
build-dyturbo-interface/leptonff_dyturbo_mass_precision \
  production.tsv migration.json 1048576 384 300000000 771 .0072973525692838
```

It uses 0.25 GeV observed bins from 50 to 116 GeV and a normalization complement
below 50. PHOTOS uses the same fixed alpha and parent domain, exponentiation
and interference, without production-dependent ME weights or pair emission.
There are no lab muon cuts. PHOTOS-off has the same parent invariant mass as
the massless analytic Born baseline. PHOTOS 3.56 settings in this campaign
have photon cutoff 1e-5 and exponentiation remainder parameter 1e-4.

PHOTOS is stratified in 0.25 GeV parent bins. Its on-minus-off control variate
preserves normalization for every sample. It also stores five polynomial
transfer moments and their covariance in each parent/observed bin, allowing
the QCD mass spectrum to be replaced without rerunning PHOTOS. A fixed parent
stratum occupies one Born bin, so region covariance can be recovered including
the negative on/off products. Independent seeds test the stated errors.

`compare_precision.py` combines those transfer moments with refined native
production and makes PNG/PDF/SVG plots, a table and a report. Supply the native
JSON, analytic JSON and independent PHOTOS JSONs. Its figure shows PHOTOS MC
errors; native input errors and fragmentation/source-grid shifts are reported
separately. It does not label these as total phenomenological uncertainties.
The local fine-bin campaign is under `results/dyturbo-mass-precision/`.

The final [numerical record](../benchmarks/dyturbo_mass_photos.json) uses 400
native parent bins of width 0.125 GeV, native relative quadrature tolerance
1e-4, 1,048,576 fragmentation nodes and 384 within-bin quadrature nodes.
PHOTOS has 2.7 billion draws across 18 independent seeds. The final reference
tightens its photon cutoff and exponentiation remainder to 1e-7 each; set
`LEPTONFF_PHOTOS_CUTOFF=1e-7 LEPTONFF_PHOTOS_EXP_EPS=1e-7` when generating it.
The 2.4-billion-draw sample at 1e-5/1e-4 is retained as a numerical control.
The plot and tables are in `results/dyturbo-mass-precision/comparison-final/`.

In 90–92 GeV, the QED-off, order-alpha, photonic-NLL and PHOTOS cross sections
are 810.553, 697.265, 703.377 and 702.570 pb, respectively. The NLL correction
is -13.223%, versus -13.322% from PHOTOS, differing by 0.0995 percentage point.
This measures a tool/model difference; it does not establish its physical
origin or certify a full fiducial EW prediction.

The largest PHOTOS MC error in an on/off ratio is 0.0937 percentage point;
in the peak bins it is below 0.0070 percentage point. The largest propagated
native-input error estimate is 0.0071 percentage point. The FFT/quadrature
and conservative source-grid refinements change ratios by at most 0.0037
and 0.0026 percentage point. The loose-to-tight PHOTOS change in the peak
rate is 0.0013 +/- 0.0160 pb (MC), consistent with zero; the per-bin difference
pull RMS is 1.05. These numerical checks do not replace scale/PDF variations,
other QED prescriptions or a wider parent production domain.

The generation manifest retains client/source/dependency hashes and a local
source/binary snapshot. The comparison can be reproduced with:

```sh
python3 adapters/dyturbo/compare_precision.py \
  --production results/dyturbo-mass-precision/native-refined.json \
  --analytic results/dyturbo-mass-precision/refined-1M-384.json \
  --photos results/dyturbo-mass-precision/tight-photos/seed*/result.json \
  --control results/dyturbo-mass-precision/refined-524k-192.json \
    results/dyturbo-mass-precision/refined-coarsened-1M-384.json \
  --photos-control results/dyturbo-mass-precision/comparison \
  --campaign results/dyturbo-mass-precision/tight-photos/campaign.json \
  --output results/dyturbo-mass-precision/reproduced-comparison
```
