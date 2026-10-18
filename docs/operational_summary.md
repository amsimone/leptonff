# Operational summary

Version 0.25 packages the generator-independent decay/tensor response as an
installable C++ library and C/Fortran interface, with configurable bins,
grouped production archives and checksum-verified multi-model replay.
MCFM now uses the common spectrum evaluator and can export production
without evaluating QED. See the [package and replay guide](packaging.md).
DYTurbo 1.4.2 native coupling is pending its source archive; the generic
spin/angular interface is ready, but no NNLO+NNLL result is claimed.

Version 0.24 adds a separate **fixed-alpha photonic NLL** mode, with explicit
unitary conditioning of the published time-like photonic kernel and one
subleading insertion globally. W+/W-/Z adapter and comparison routing are
available as `--photonic-nll` / model `photonic_nll`. See the
[normalisation, matching and validation contract](photonic_nll.md).
The frozen v0.23 high-level comparison plots still contain photonic LL;
this increment does not launch another production campaign.

The NLO QCD W±/Z fragmentation benchmark is operational. It can produce
bare-muon fiducial rates, observed differential spectra, PHOTOS comparisons,
normalised bin fractions and full grouped-event covariance. The coupled
one-family reference, variable-flavour leptonic NLL calculation, sourced-HVP
pointlike hadronic benchmark with LO QCD, and fixed-alpha photonic comparison
are supported modes. This completes the NLO W/Z benchmark and
differential-measurement milestone. The intended accuracy is NLO QCD
production with final-state NLL QED in the specified flavour/input model.
Mixed QCD–QED splitting kernels and fitted hadronic fragmentation are optional
extensions, not prerequisites to that benchmark. Numerical precision and
alignment to a selected measurement still need to be demonstrated; NNLO QCD
remains a later programme. See the [accuracy audit](qed_accuracy_audit.md).

## What is implemented

- Distribution-aware O(alpha) boundaries, complete one-family time-like P0/P1,
  coupled NLL evolution, running coupling and explicit expansion coefficients.
- Full charged-flavour time-like P0/P1 with separate produced-colour and
  closed-loop factors, full-parent evolution of selected observed columns,
  finite heavy-flavour and photon matching, two-loop running and leptonic
  on-shell/MSbar input conversion. The e/mu/tau preset includes the tau threshold.
- Spin-correlated, unitarily normalised finite decay matching for W±, including
  charged-W photon emission, and the charge-symmetric neutral-current sector.
- Muon and photon hard-parent channels in the coupled match, actual
  lattice-overlap subtraction, continuum NLO restoration and threshold-consistent
  hard-scale coupling and second-order expansion.
- Conservative signed resolved-fraction cells, exact endpoint atoms and exact
  cell integrals for two-leg mass/pT cuts. Atomic references remain supported.
- An isolated MCFM 10.3.1 NLO QCD adapter preserving signed real/dipole grouping,
  original spin/recoil information and shared radiation draws on counterevents.
- PHOTOS 3.56 decays on the same production contributions, with separate
  massive Born mapping, conservation checks and deterministic common seeds.
- Observed W muon pT, absolute eta and mT spectra; observed Z pair mass,
  minimum muon pT and maximum absolute muon eta spectra. Overflows are retained.
- Fixed-alpha photonic LL plus one finite insertion for the coupling/pair-aligned
  PHOTOS benchmark. The auxiliary `ff2` column in this mode is FF-only.
- Full event covariance across 136 spectrum columns and eight total-rate
  columns, normalised shapes, paired PHOTOS differences and independent
  reconstruction of native host errors.
- Separate seven-point QCD-scale and three-point fragmentation-scale studies,
  grid/cutoff/seed diagnostics, per-run provenance and PNG/SVG plots.

The [measurement and covariance contract](lab_spectra.md) states exactly what
each spectrum counts. The [v0.20 benchmark](../benchmarks/lab_spectrum_reference.md)
records checks, numerical results, warnings and reproduction commands.

## Verified scope

The core suite and independent numerical references cover distributions,
moments, conservation laws, time-like kernels, evolution, matching expansions,
spin limits, endpoint/cell integration and laboratory measurements. Earlier
independent RECOLA validation covered 486 polarised real-current points.
The v0.21 suite adds all-flavour momentum/valence sums, charge/colour controls,
finite threshold maps, inverse coupling transport and cubic matching-remainder
checks. An independent analytic-source/exponential-midpoint calculation checks
the full-parent moment results and selected-column x-space moments.

The new spectrum tests compare against literal signed lattice tensor sums,
direct observed kinematics and the original photonic geometric operator.
The statistical test includes signed events, zero calls, bin migration,
normalisation correlations and rejection of an incomplete covariance file.
At each actual host iteration, event sums and covariance diagonals must agree
with MCFM's canonical histograms. Every spectrum's integral and covariance
must reproduce its total rate. Zero-coupling, grid, cutoff, NLO scale-cancellation
and independent LL/NLL-backend controls are also audited.

Passing these checks establishes the implementation contract. Some reduced-call
runs have high integration chi-squared, and rare-tail/angular fluctuations
remain visible. Those warnings are preserved. Scale envelopes and grid shifts
are reported with their sampling errors; none is promoted to a precision bound.

## What is still needed for phenomenology

| Gate | Current position | Required work |
| --- | --- | --- |
| Aligned tool benchmark | Completed W pT/|eta| and Z mass comparison; two higher-angular seeds for W+/Z and ten variance-control seeds for W− | Improve the stated statistical/numerical limits for a more demanding precision target; validate production recoil for precision data |
| Pure-QED NLL flavour/input prescription | e/mu/tau audited; sourced HVP/pointlike quark prescription with strong FF evolution off and measured input/threshold sensitivities | State the chosen infrared/pair prescription for data; resolve numerical grid sensitivity before precision flavour comparisons |
| Hadronic vacuum polarisation | Sourced PDG input, alternative installed Jegerlehner evaluation and explicit fermionic MSbar matching | Higher conversion/mass terms and input correlations need the chosen precision treatment; HVP alone does not determine real fragmentation |
| Experimental lepton/pair definition | Inclusive tagged W leptons and distinct-parent opposite-charge Z pairs | Implement the experiment's extra-pair, leading-pair, charge, isolation and prompt-muon conventions where they differ |
| Differential Z scope | Charge-symmetric pair mass and extrema | Add the requested charge-summed single-muon spectra or full parity-odd spin information for charge-asymmetric observables |
| Numerical/statistical precision | Full covariance, seed ensembles and grid/domain/cutoff controls; fixed-resolution domain changes negligible | Grid responses reach 0.6% of Born bins in the coarser HVP Z control; geometry/cutoff sampling remains noisy. No per-mille certification |
| Production uncertainties | Seven QCD scale points and fragmentation-scale study available | Use the chosen PDF order/input set and evaluate PDF/alpha_s uncertainties; validate recoil resummation near Jacobian/low-recoil regions |
| Data alignment | Published bare/post-FSR candidates identified | Match collision energy, fiducial cuts, line shape, widths/branching conventions and published bin definitions; the current 13 TeV setup is not the 8 TeV CMS measurement |
| Further accuracy | NLO QCD production and factorised final-state QED | NNLO QCD host integration, full NLO electroweak effects and exact mixed corrections remain later milestones |

The [new comparison](../results/qed_comparison_v023/comparison.md) passes 144
main checks with no host warning flags. Its 15 auxiliary HVP-input checks
retain five reduced-call W+ flags. Eleven focused tests and 48 default
adapter-equivalence checks pass. The [benchmark](../benchmarks/qed_comparison_reference.md)
records achieved errors, numerical limits and the variance-control study.

The v0.22 increment passes 60 core tests and 101 host checks over 39 control
runs. Thirteen reduced-call W+ convergence flags remain in the audit.
Independent source/analytic-moment/exponential checks and a read-only installed
HVP routine validate the new inputs and evolution. These are implementation
checks, not differential precision certification. Mixed kernels and fitted
hadronic boundaries concern the additional combined evolution programme;
they do not block the pure-QED comparison.

The [v0.21 flavour/threshold contract](flavour_threshold_matching.md) and
[benchmark](../benchmarks/flavour_threshold_reference.md) document the new
mode, its independent moment/evolution checks and W±/Z controls. Full flavour
mixing is implemented, including feedback through other charged-parent rows;
the selected legacy 3×3 output is not a closed momentum system. Generic custom coloured entries retain the pure-QED reference. The new
[hadronic contract](hadronic_qcd_matching.md) adds LO QCD and a sourced HVP
anchor in a separately named pointlike-boundary benchmark. It preserves
matching-scale and input provenance without treating that boundary as fitted.

For the coupling/pair-aligned analytic-versus-PHOTOS comparison, use the
photonic branch. For leptonic NLL running and pair effects, use
`--coupled-hard --fraction-cells --flavours leptons`; the default numeric alpha
is then treated as on-shell input and converted at mu0. The latter comparison
to PHOTOS is diagnostic. A quantitative precision-data claim requires the open
gates above.
