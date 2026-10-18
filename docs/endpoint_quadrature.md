# Endpoint geometry and measurement-aligned FF panels

Version 0.14 improves the numerical action of the existing one-family matched
single-muon spectra. It introduces no new hard coefficient or logarithmic
accuracy claim. The [spectrum measurement](lepton_spectra.md), inclusive pair
counting, common MSbar coupling and physical limitations are unchanged.

## Origin of the panel sensitivity

A muon from a parent with energy at most M/2 cannot pass an extra pT threshold
t when its FF fraction satisfies `f <= 2*t/M`. The finite geometry therefore
vanishes when both observed fractions lie below this boundary. Interpolating
across that contact can mix a forbidden region with an allowed one. Agreement
of cubic FF moments and a small total-rate shift do not diagnose this error.

The largest v0.13 changes were near the 45 GeV edge of the pT spectrum. A
second source was compression of the endpoint/resolved geometry, which enters
the hard dressing already at O(alpha²). The new calculation addresses both.

## Direct endpoint strips

Write the weighted FF source arrays as `a_j = z_j D_a(z_j) dz` and similarly
for b, with `z_0=1`. Let `q_a` and `q_b` be the existing linear signed cubic
projection onto fractions `f_i`, keeping the lattice endpoint in index zero.
For a finite geometry G, the new action is

```text
a_0 b_0 G(1,1)
+ b_0 sum_{j>0} a_j G(z_j,1)/z_j
+ a_0 sum_{j>0} b_j G(1,z_j)/z_j
+ sum_{i,k>0} q_a[i] G(f_i,f_k) q_b[k].
```

These four disjoint sectors partition the numerical source. Both endpoint/
resolved strips are evaluated directly on the original FF lattice; only the
resolved/resolved sector uses compression. Muon geometry is symmetric, while
the two photon strip orientations are stored and applied separately. Endpoint
weights may have either sign. No positivity clipping, reweighting or bin
normalisation is introduced.

The same action is used in the resummed terms and their explicit expansion.
The geometry in C1 uses G(1,1). The geometry in C2 involves one first-order FF
and one Born endpoint, so it now has no FF-panel interpolation error. The
compressed interior first contributes to the hard-matched rate at O(alpha³).
The reported C2 remains a finite-grid coefficient; removing panel dependence
does not remove its grid or hard-integration errors.

`CoupledDecayOptions.exact_endpoint_strips` enables this action in the general
fiducial operator. `CoupledLeptonSpectrum` always enables it, as it already
enables conservative mass-cell integration. The other reference drivers keep
their previous numerical prescription by default. Direct strip evaluation
adds work linear in the retained FF prefix, rather than evaluating the full
quadratic lattice tensor. "Direct" still includes the stated real-integration
tolerances and omitted soft/collinear slivers.

## A common panel partition

`CoupledDecayOptions.ff_fraction_breaks` supplies additional fraction-space
panel boundaries. The uniform square-root-loss partition remains the baseline;
new boundaries are inserted into that partition. Each active subpanel uses
four fixed interpolation nodes. Empty subpanels are removed from the source
projection. The endpoint stays explicit and cubic moment preservation remains
valid for signed weights.

For simultaneous pT and eta spectra use

```cpp
auto options = lepton_spectrum_options(M, cuts, pt_edges, eta_edges, base);
CoupledLeptonSpectrum pt(M, cuts, LeptonObservable::transverse_momentum,
                         pt_edges, grid, options);
CoupledLeptonSpectrum eta(M, cuts, LeptonObservable::pseudorapidity,
                          eta_edges, grid, options);
```

The helper inserts every physical `2*t/M` pT contact and the eta/pT angular
transition `2*pt_min*cosh(eta_edge)/M`. It selects direct endpoint strips and
mass-cell integration too. Both histograms and all their cumulatives use the
same partition, support and source projections. Refining the baseline panel
count retains all physical contacts. The reference driver uses this helper.
For custom separate calculations, use the same prepared options wherever a
common finite-grid pair total is required.

These boundaries refine the numerical integration. They neither move the
original FF lattice nor smear pT/eta measurements. In particular, each
original pair cut and the two-entries-per-pair contract stay in force.

## Audit and reproduction

`--dump` now includes the nonuniform panel edges, interior geometry matrices
and direct endpoint arrays. Endpoint arrays contain `G(z_j,1)/z_j` or
`G(1,z_j)/z_j`, matching the weighted source convention. Each production point
reports its active FF node count, and the header records the common contacts
and geometry prescription. The independent script rejects mismatched controls
or panel-dependent C2 values in a panel refinement. CTest writes a checked
coarse artifact into the build directory. Optional `--checked-coarse` reuses
that audit only when source, executable and both configuration hashes match;
its hash and original audit provenance are retained in the final manifest.

```sh
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 2 --scan > results/lepton_spectrum_endpoint_scan.json
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 8 > results/lepton_spectrum_endpoint_panels8.json
python3 python/check_lepton_spectrum_reference.py build/leptonff_lepton_spectrum \
  configs/reference.conf configs/lepton_bins.conf \
  --checked-coarse build/lepton_spectrum_endpoint_checked_coarse.json \
  --input results/lepton_spectrum_endpoint_scan.json \
  --panel-input results/lepton_spectrum_endpoint_panels8.json \
  --output results/lepton_spectrum_endpoint_reference.json
```

The C++ support test uses a coarse lattice for which only the endpoint can
produce a 45 GeV muon. It compares full signed muon/photon geometry against a
literal lattice tensor, checks panel independence and verifies the cubic
coupling remainder after C2 subtraction. Existing bin closure, scale, mass-log
and massive-reference checks still apply.

Python independently reconstructs the source arrays, the nonuniform cubic
projection, both endpoint strips by allowed-azimuth integration, the mixed
geometry algebra and C1/C2. It checks that a forbidden lattice fraction never
projects onto an allowed pT region. Interior matrices are sampled by separate
angular quadrature and used for the remaining contraction audit; this is not
an independent angular reconstruction of every matrix entry. Numerical changes
and the remaining grid sensitivity are recorded in the
[benchmark](../benchmarks/endpoint_quadrature_reference.md).
