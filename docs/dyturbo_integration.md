# DYTurbo 1.4.2 integration status

Target source:
`https://dyturbo.hepforge.org/downloads/dyturbo-1.4.2.tar.gz`.

The common library, direct production-tensor response, explicit angular-weight
conversion and recoil-frame helper are implemented and tested. **There is
not yet a source-verified native DYTurbo adapter or an NNLO+NNLL production
result.** No existing plot should be relabelled as such.

DYTurbo is not vendored. Obtain the requested release separately, then prepare
an isolated source tree:

```sh
wget https://dyturbo.hepforge.org/downloads/dyturbo-1.4.2.tar.gz
python3 adapters/dyturbo/prepare_source.py dyturbo-1.4.2.tar.gz
```

This creates an isolated source tree under
`external/dyturbo-1.4.2-work`, pins the archive SHA256, confirms the release
string and lists real source candidates for the angular/integrand interfaces.
It refuses overwrites, traversal and escaping links. Its automated tests
use synthetic tar fixtures; they are not DYTurbo integration tests.

## Native coupling that remains to implement

The published DYTurbo method factorises production from leptonic angular
acceptance and uses fast angular quadrature/analytic integrals. The radiative
response needs to enter those acceptance integrals, consistently in the
resummed, asymptotic/counterterm and finite-order pieces. Applying a scalar
FSR factor to completed one-dimensional histograms is not this coupling.
See the [DYTurbo paper](https://arxiv.org/abs/1910.07049).

The source-level implementation must:

1. Identify the actual 1.4.2 production/acceptance interface, normalisation,
   recoil prescription and angular axes. Map them explicitly to the common
   decay tensors; retain signed/zero production terms without Born division.
2. Establish QED-off closure with native DYTurbo, first at matching NLO
   settings and compatible W-width/EW/branching inputs, then at NNLO+NNLL.
3. Apply one radiative measurement to every matched QCD contribution and
   preserve its cancellation and cutoff controls. The current W decay is
   a factorised/pole-type correction; it is not a full off-shell EW amplitude.
4. Cache production coefficients or nodes with PDF, QCD and resummation-scale
   variations. Replace expensive repeated angular acceptance by validated
   decay responses where possible. Interpolation and quadrature errors need
   their own convergence study, separate from photon Monte Carlo variance.
5. Validate photonic order-alpha and LL/NLL limits, then compare the requested
   bare-lepton W pT/eta and Z mass spectra with a shared PHOTOS reference.

The result would be NNLO+NNLL QCD combined with the declared matched NLL QED
decay prescription. It would not supply the complete mixed QCD-EW correction.
Jets or QCD-radiation-sensitive observables are outside DYTurbo's inclusive
production interface and the present decay-only replay schema.
