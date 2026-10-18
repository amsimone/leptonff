# Matched single-muon pT and signed eta spectra

Version 0.14 refines endpoint geometry and inserts measurement contacts into a
common FF panel partition; see the [numerical contract](endpoint_quadrature.md)
and [current benchmark](../benchmarks/endpoint_quadrature_reference.md).

Version 0.13 adds `CoupledLeptonSpectrum` and `leptonff_lepton_spectrum` to the
[coupled fiducial reference](coupled_fiducial_matching.md). Both muons fill the
histogram for every selected opposite-charge pair. Every additive component,
including C1 and C2, therefore sums to **twice the pair rate**. The `total` field
is the pair rate, `sum` is the histogram sum, and closure is `sum - 2*total`.
These integrals count identified muons in inclusive distinct-parent pairs;
they are not exclusive two-muon event probabilities.

The original pair selection stays in place: both muons pass the fiducial pT
and absolute-eta cuts, and their pair mass lies in the original mass window.
An extra cumulative threshold acts on one muon at a time. In particular,
raising both pair pT cuts would incorrectly remove pairs in which only one
muon exceeds the histogram threshold.

## Measurement and matching

Let F be the original pair measurement. The mean one-lepton cumulatives are

```text
P(t) = F * [ I(pT_1 > t) + I(pT_2 > t) ] / 2
A(t) = F * [ I(|eta_1| < t) + I(|eta_2| < t) ] / 2.
```

For pT bin `(lo, hi]`, the two-entry integral is `2*[P(lo)-P(hi)]`.
pT edges start at `pt_min` and cover the kinematic endpoint M/2. Thresholds
at or above M/2 give zero. At fixed parent mass there is no physical DY tail
above that endpoint; a bin crossing it contains the integrable Born Jacobian
peak. pT is continuously averaged over rotations, so there is no pT endpoint
atom to assign to a particular bin.

Signed eta uses global reflection symmetry of the unpolarised parent at rest.
A positive bin has integral `A(hi)-A(lo)`, a negative bin
`A(-lo)-A(-hi)`, and a bin crossing zero `A(-lo)+A(hi)`. Edges must partition
`[-eta_max,+eta_max]`; asymmetric interior edges and bins crossing zero are
allowed. This shortcut must be replaced for a boosted or spin-correlated host.

Each cumulative retains the same FF support, hard source arrays and linear
signed cubic FF projection as the pair operator. Endpoint/resolved finite
geometry is evaluated directly on the original FF lattice; the cubic
projection is used only when both fractions are resolved. The selected angular action
is the average of two rotational band intersections, restricting only the
chosen muon's band each time. Both muon-parent orientations and photon-parent
channels receive the same measurement. Projected hard terms and the first/
second expansion use that measurement too. The complete continuum O(alpha)
target is restored with the existing additive overlap prescription.

The Born pT lattice action depends on both fractions; a truncated symmetric
literal pair sum evaluates it. The Born eta action retains the fast prefix
sum because the two directions are back-to-back. The conservative **mass**
cell rule from the [mass spectrum](mass_spectrum.md) is retained for both
histograms. There is no additional pT or eta cell smearing. Finite geometry is
continuously integrated at the shared FF quadrature fractions. Signed weights
are preserved, with no clipping or bin renormalisation.

## Output and variations

Bin integrals are dimensionless relative to the massless inclusive vector Born
normalisation. pT densities divide by bin width in GeV; eta densities divide by
its dimensionless width. Couplings and compression diagnostics are shared
metadata. Each point reports both spectra and a common pair total, including
FF, muon/photon hard terms, raw and matched rates, fixed order, and C1/C2.

The model remains one-family QED with an unpolarised pure-vector parent at rest,
fixed M, and inclusive identified opposite-charge pairs across distinct hard
parents. It omits physical flavours/thresholds, axial and spin effects, boosts,
production PDFs and an NNLO host adapter. It makes no physical precision claim.

`--scan` doubles the FF panels, doubles the grid twice, doubles the domain at
fixed finest spacing, and divides hard slivers/tolerances by four. It then
varies fragmentation and starting scales by factors 1/2 and 2. Each variation
acts on both spectra together; hard alpha remains at M. Starting-scale changes
transport the common MSbar input before recomputing evolution and boundaries.
Report these correlated shifts separately from numerical refinements. Unknown
constant two-loop matching terms and the missing physical content are not
estimated by the narrow total-rate band.

```sh
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 2 --scan > results/lepton_spectrum_endpoint_scan.json
./build/leptonff_lepton_spectrum configs/reference.conf configs/lepton_bins.conf \
  8192 8 > results/lepton_spectrum_endpoint_panels8.json
python3 python/check_lepton_spectrum_reference.py build/leptonff_lepton_spectrum \
  configs/reference.conf configs/lepton_bins.conf \
  --input results/lepton_spectrum_endpoint_scan.json \
  --panel-input results/lepton_spectrum_endpoint_panels8.json \
  --output results/lepton_spectrum_endpoint_reference.json
```

The bin file requires unique `pt_edges` and `eta_edges` comma-separated entries.
`--dump` additionally exports the retained weighted FFs, hard projections,
shared fractions and selected hard-geometry samples for an independent audit.
Saved checks hash source, executable, input and both configuration files.

## Verification scope

C++ tests check Born bin integrals and the Jacobian endpoint, a one-leg threshold
with only one qualifying muon, a signed non-endpoint literal tensor, two-entry
closure, eta reflection, scale cancellation in C1 and absence of retained
alpha² mass-log dependence in each C2 scale shift. Selected pT and eta
cumulatives also approach the exact massive real/virtual reference.

Python independently builds FFs by causal series and hard source cells by
z quadrature, then evaluates literal signed lattice pairs. It checks the FF
and projected hard bin terms. Direct allowed-azimuth quadrature checks selected
muon/photon geometry samples; the leading-power massive spectrum plus direct
geometry supplies a separate continuum O(alpha) bin reference. This samples,
rather than reconstructs, the full resummed geometry matrices. Inherited
coupled tensor checks and new bin-specific grid/panel/sliver scans supply the
remaining numerical evidence. See the [benchmark](../benchmarks/lepton_spectrum_reference.md).
