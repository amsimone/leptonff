# Laboratory FF convolution and NLO QCD host pilot

This is the retained v0.15 FF-only contract. The optional v0.16 finite decay
match is documented [separately](spin_decay_matching.md).

Version 0.15 adds a W one-leg and Z two-leg laboratory measurement, a working
MCFM 10.3.1 analysis adapter, and a common-production PHOTOS comparison.
This increment is **Born-parent FF-only convolution on NLO QCD production**.
It does not supply the finite W/Z O(alpha) hard-decay coefficient and therefore
does not establish a fully matched physical NLL cross section.

## Measurement and overlap contract

`LabFragmentation` consumes host massless lepton four-momenta, including their
spin-correlated angular distribution and QCD recoil. It never regenerates an
isotropic decay. The W action rescales the charged momentum by z, leaving the
neutrino unchanged; the visible collinear remnant is not missing momentum.
W pT, eta, neutrino pT and transverse-mass cuts are applied to the resulting
measurement. At fixed input momenta, `mT^2(z) = z*mT^2(1)`.

The W pilot observes the charge of the primary hard lepton. It is a tagged
W-charge sample, not an exclusive one-muon event definition; migrating
wrong-sign leptons from pairs between observed-charge samples remains a
physical-content extension. The Z action includes same/same and opposite/
opposite identified-pair terms from the two hard directions, as in the existing
reference contract. No extra-pair veto is imposed.

`SignedFractionTable` forms signed number-weight prefixes from the original
energy-weighted FF lattice, dividing by z only on supported measurements.
The unresolved endpoint atom is retained, bounds are strict, and no negative
weights are clipped or normalised. Queries below the explicit positive support
limit fail. The pilot uses z_min=1e-4, safely below the support required by
pT>25 GeV at 13 TeV. Resolved bins remain lattice atoms; grid refinement is
required for sharp measurements.

For each host contribution, the reported result is

```text
FF_matched = FF_resummed + a_initial*(C1_FF_continuum - C1_FF_lattice)
```

The continuum whole-plus action is integrated analytically at the laboratory
cut bounds. Its implementation is independently checked against distribution
quadrature. The actual lattice overlap is removed, not an unrelated continuum
approximation to it. C2 is also extracted, including the two-leg C1 product.
This matches the **FF contribution** through O(alpha); a process-dependent
finite hard correction is still needed. A QCD host alone does not supply it.

## MCFM integration

The adapter replaces the old-style `nplotter_user` object in an isolated
executable, reusing compatible local GCC objects/libraries without modifying
the host installation. `build_adapter.py` records source/library/executable
hashes. The dependency-free core remains separately buildable and testable.

Use `histogram.newstyle=false`: the new-style dispatcher bypasses this hook.
The runner verifies that the FF callback actually ran and produced its required
histograms. MCFM configuration keys must be indented beneath their sections;
the runner writes that syntax explicitly. Native lepton cuts are disabled for
the convolution, so migration is applied before observed cuts. Off-shell
Z/gamma* parents extend above the observed 120 GeV mass boundary; parents below
60 GeV cannot migrate upwards under the collinear FF action.

Each convolved observable is evaluated before calling `bookplot`. The host
temporarily sums real and dipole contributions per integration point before
squaring them. The adapter advances the histogram index on every callback.
Small mass-shell roundoff residuals from mapped momenta are repaired in energy
only, preserving three-momenta; residuals above a relative energy tolerance of
1e-6 fail, and the actual maxima are printed. This is numerical repair of a
declared massless host, not a massive-lepton approximation.

Seven rate components and eta histograms are recorded. The Z eta histogram
counts the first parent-direction entry per selected pair; it is not the
existing both-muon histogram. Dedicated observed-pT/mass histograms in this
host adapter remain to be implemented.

The host combines iterations independently for each histogram/bin. Consequently
printed means need not obey algebraic closure exactly. Directly accumulated
`ff_shift`, `ff_beyond_order1` and `photos_minus_ff` histograms provide paired
errors; subtracting printed central rates does not reproduce those estimators.
Full cross-bin covariance and a joint uncertainty prescription remain open.

## PHOTOS comparison

The optional adapter uses the local PHOTOS 3.56 HEPEVT build. It is restricted
to one OpenMP thread because PHOTOS has global state. Both branches receive
the same signed NLO QCD contributions, before observed cuts. Real and mapped
dipole contributions use a photon random stream seeded from the original real
lepton momenta; independent counterevent photon draws would spoil cancellation.

PHOTOS runs as a standalone boson decay with exponentiation and interference
on, W matrix-element correction on, Z matrix-element correction off, alpha(0)
=1/137.035999084 and photon cutoff 1e-5. No QCD shower or second FSR source is
attached. Massless host decays are mapped to physical-muon-mass two-body decays
in the parent rest frame, preserving the parent and decay direction, then
boosted back to the lab. The separate `photos_born` measurement exposes this
mass projection. Momentum conservation and common-seed repeatability are
tested for both W charges and Z.

This is an explicit generator configuration, not a universal PHOTOS accuracy
claim. The current FF branch has a one-family MSbar coupling/boundary and pair
content; PHOTOS supplies photonic decay radiation with its stated finite terms.
There is no finite conversion or hard-decay alignment yet. Differences between
the branches therefore include these convention/content differences and must
not be attributed to NLL resummation alone.

See the [benchmark and reproduction commands](../benchmarks/nlo_host_reference.md)
for the actual runs. The next physics task is the W-specific spin-correlated
finite hard match, followed by the corresponding Z match, aligned photonic/
flavour/coupling inputs and expanded-order generator comparisons. NLO W/Z
integration and tool comparisons precede the NNLO QCD extension.
