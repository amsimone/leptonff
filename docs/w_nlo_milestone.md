# First hadronic milestone: bare-muon W production at NLO QCD

Version 0.17 now supplies photonic LL dressing of the W±/Z hard operator and
boundary with one finite insertion globally; see the
[dressing contract](photonic_lab_dressing.md) and
[benchmark](../benchmarks/photonic_dressing_reference.md). Coupled physical
NLL and precision phenomenology remain subsequent milestones.

17 October 2026. Version 0.15 now implements the NLO W±/Z host FF-only pilot and
common-production PHOTOS comparison. See the [implemented contract](lab_host_matching.md)
and [benchmark](../benchmarks/nlo_host_reference.md). The finite W hard-decay
match and physical-content alignment were the next steps at that version.
Version 0.16 now provides additive spin-correlated W matching, a massive W
anchor and charge-symmetric Z matching, with a separate fixed-alpha photonic
NLO+LL comparison branch. See [contract](spin_decay_matching.md) and
[benchmark](../benchmarks/spin_decay_reference.md). Full NLL hard dressing,
physical multi-flavour content and the phenomenology error budget remain open.
W2 is now implemented at additive O(alpha); W3 has a restricted photonic
comparison, with numerical component precision and broader physical content
still under study.

## Deliverable and order of work

Produce separate pp → W+ → μ+ν and pp → W− → μ−νbar fiducial cross sections,
muon eta distributions and charge asymmetry, then coarse muon pT and transverse
mass distributions. Use fixed-order NLO QCD production as the first host.
The existing FF kernels, signed distribution action, evolution, moment and
expansion checks are reusable. The rest-frame two-muon geometry is not.

| Step | Deliverable | Acceptance |
| --- | --- | --- |
| W0 | Native NLO host and identity measurement | Native rates and bins reproduced with correlated integration errors; inputs, units and generation cuts recorded |
| W1 | One-leg expanded O(α) and resummed LL FF pilot | Signed endpoint/plus action; alpha expansion, cut migration and subtraction-parameter checks; common quadrature across correlated terms |
| W2 | W-specific O(α) hard decay plus FF match | Spin-correlated real/virtual/subtracted decay, finite terms, independent massive/fixed-order comparison, overlap removal and scale cancellation |
| W3 | Physical-content and comparison benchmark | Documented flavours/coupling/width convention; aligned reference; separate numerical, QCD/PDF and fragmentation variations |

W1 is useful integration evidence before the full physical flavour extension.
It must be labelled NLO QCD plus FF-only collinear radiation in its stated
theory. A matched NLL cross-section claim requires W2 and the W3 theory inputs,
including the required hard-parent channels. The current one-family kernels
alone do not provide a physical full-QED result.

## Laboratory measurement

Use one charged fragmentation leg, with observed momentum `p_mu = z*k_mu`.
The neutrino and QCD momenta are unchanged in this collinear FF action. For the
initial truth definition, missing transverse momentum is the neutrino vector.
The collinear visible remnant `(1-z)*k_mu` must not be reassigned to missing
momentum. Explicit hard-photon recoil belongs to the W2 real-decay mapping.

Define

```text
mT^2 = 2*pT_mu*pT_nu*(1-cos(phi_mu-phi_nu))
pT_mu(z) = z*pT_mu(1), eta_mu(z) = eta_mu(1)
mT^2(z) = z*mT^2(1)                 [pure collinear FF action only]
```

Draft engineering setup: 13 TeV pp, pT_mu > 25 GeV, |eta_mu| < 2.4,
pT_nu > 25 GeV and mT > 40 GeV, with no photon recombination, isolation or
extra-lepton veto. These are proposed validation cuts, not an adopted data
comparison. Freeze beam energy, bins, PDFs/alpha_s, CKM inputs, masses/widths,
EW scheme and decay normalisation before producing a comparison table.
Use charge-labelled inclusive identified-muon counting if pair channels are
enabled; an exclusive one-muon event definition needs separate justification.

In the pilot, retain the host's leptonic spin correlations and QCD recoil.
Apply cuts in the convolved measurement. Move or relax native cuts as needed
and demonstrate generation-cut independence. Validate covariance across FF
nodes, real/counterevents and scale variations; individual squared weights
are not a substitute for correlated integration errors.

The W hard decay needs its own charged-current calculation. The
[pole-approximation treatment](https://arxiv.org/abs/1511.08016) explains the
gauge-consistent separation with radiation from the charged W. Reusing half
the neutral-current hard coefficient does not supply this matching.

## First host probe: local MCFM 10.3.1

Read-only inspection found the following under
`/Users/simone.amoroso/Work/generators/MCFM-10.3.1`:

- A built `build-gcc/mcfm` executable and an earlier LO smoke-test artifact.
  Neither establishes a successful NLO W run for this project.
- W input examples and `src/User/nplotter_W_only.f`. The routine receives
  `(px,py,pz,E)` momenta and signed weights; its switch distinguishes real
  contributions from subtraction terms. W+ uses charged-lepton index 4;
  W− uses index 3. Check the neutrino assignment and flavour convention in
  the selected process before constructing records.
- `src/User/nplotter_user.f90` receives the dipole number `nd`;
  `src/Procdep/realint.f` and `virtint.f` call plotting with mapped momenta.
  These are candidate insertion points, not a verified correlated API.
- `nplotter.f` fills its total-cross-section histogram before the user hook.
  A new convolved fiducial rate therefore needs its own histogram/accounting.
  Native selection and integration-error accumulation still need tracing.
- The generic W+ input uses `part=resNLO`; an NLO study input uses
  `removebr=true`. Neither is a ready leptonic fixed-order baseline. Prepare
  a dedicated `part=nlo`, `removebr=false`, `ewcorr=none` configuration,
  consistent with the chosen width treatment and accessible PDF set.

The [official W documentation](https://mcfm.fnal.gov/processes/process1.html)
lists processes 1 and 6 and NLO dipole subtraction. First exercise the existing
installation without modifying it. Develop any source adapter in an isolated
copy within this workspace, with a pinned source digest and reproducible patch.
Core tests must continue to run without the host dependency.

## Comparison contract

The primary comparison is the bare-muon FSR correction and fiducial rates from
PHOTOS/PYTHIA versus analytic fragmentation on common NLO QCD production. This
precedes experimental-data comparisons. The same host/measurement machinery
must also pass NLO QCD Z-cross-section tests before an NNLO extension.
Independent fixed-order and photon-resummation references need
aligned inputs, perturbative content and radiation settings; do not duplicate
FSR. Compare expanded O(α), LL and matched NLL only with their hard terms stated.

Before selecting experimental tables, confirm whether muons are bare post-FSR,
dressed or unfolded to pre-FSR Born kinematics. For example, the
[ATLAS 8 TeV W cross-section measurement, section 3.1](https://arxiv.org/pdf/1904.05631)
defines its fiducial quantities at Born level. It cannot directly validate a
bare-muon prediction. The CMS audit below identifies both an integrated
post-FSR muon-only W benchmark and a bare-level W shape. Validate the exact
analysis convention in the host setup before adopting either as a comparison.

Prioritise rates and eta bins over narrow Jacobian regions. NLO QCD is an
integration and physics benchmark; a precision W-mass fit needs a separate
assessment of recoil/resummation and other corrections. Full NLO EW, exact
mixed corrections and NNLO QCD remain later deliverables, with explicit claims.

## CMS bare-muon data audit (17 October 2026)

- **Confirmed W shape:** [CMS SMP-14-012, 8 TeV](https://cms-results.web.cern.ch/cms-results/public-results/publications/SMP-14-012/),
  JHEP 02 (2017) 096, arXiv:1606.05864. Additional Table 1 specifies a
  normalized post-FSR W transverse-momentum spectrum with zero photon dressing.
  [HEPData version 2](https://www.hepdata.net/record/ins1471281) includes it as
  Supplementary Table 1, DOI `10.17182/hepdata.77221.v2/t8`.
  The main tables are pre-FSR. The muon fiducial thresholds in section 5 are
  pT > 20 GeV and |eta| < 2.1. Check auxiliary qualifiers and its post-FSR
  momentum definition before coding a comparison: this is not a muon-pT spectrum
  or an absolute normalization, and the low transverse-momentum region needs
  separate QCD recoil/resummation control.
- **Integrated muon-only W benchmark:** [CMS SMP-12-011, 8 TeV](https://arxiv.org/pdf/1402.0923),
  appendix A table 5, gives W+ → μ+ν, W− → μ−νbar and their sum as
  `3.10 ± 0.01 ± 0.04 ± 0.08 nb`, `2.24 ± 0.01 ± 0.02 ± 0.06 nb`, and
  `5.34 ± 0.02 ± 0.06 ± 0.14 nb`, respectively; errors are statistical,
  systematic and luminosity. Muon cuts are pT > 25 GeV and |eta| < 2.1.
  Section 8 of [arXiv:1606.05864](https://arxiv.org/pdf/1606.05864) explicitly
  identifies the earlier measurement's fiducial convention as post-FSR.
  Its quoted 5.42 nb earlier result is the electron/muon combination, not
  the muon-only 5.34 nb. [HEPData](https://www.hepdata.net/record/ins1280200)
  tables 1–3 contain charge-separated and summed absolute W results.
  Use the fiducial muon entries; keep full-phase-space and combined-channel
  numbers separate. The paper uses missing momentum for fitting yields and
  imposes no minimum reconstructed missing-transverse-energy requirement.
- **Confirmed absolute dimuon benchmark:** [CMS SMP-13-003](https://arxiv.org/pdf/1310.7291),
  JHEP 12 (2013) 030, table 6, reports a 7 TeV post-FSR detector-acceptance
  Z/gamma* → μμ cross section of
  `516.5 ± 0.4 (stat) ± 4.9 (exp syst) ± 1.1 (theory syst) ± 11.4 (lumi) pb`
  for 60 < m_mumu < 120 GeV. Section 4.1.1 specifies leading/trailing pT
  thresholds of 14/9 GeV, |eta| < 2.4 and no photon recovery. Table 7 and
  section 5.2 also supply normalized post-FSR mass and mass/rapidity results.
  Preserve the analysis pair/normalization definitions when using these data.
- **Not bare W benchmarks:** [CMS 5.02/13 TeV inclusive results](https://arxiv.org/pdf/2408.03744),
  section 5, and [CMS 13.6 TeV inclusive results](https://arxiv.org/pdf/2503.09742),
  section 4, define generator fiducial leptons before FSR. The
  [13 TeV W double-differential measurement](https://arxiv.org/pdf/2008.04174)
  uses dressed leptons for those fiducial spectra. Detector muon tracks alone
  do not establish the published particle-level convention.

The v0.18 [coupled hard increment](coupled_lab_matching.md) now connects the
physical W±/charge-symmetric Z spin operators to one-family NLL evolution
and all hard-parent channels, with common-coupling and C2 checks. Physical
flavours/thresholds and a precision phenomenology budget remain open.
