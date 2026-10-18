# Evolution backend decision — reference stage

Decision for version 0.2.0: use the standalone fixed-coupling non-singlet
log-loss solver described in [ll_evolution.md](ll_evolution.md). The implementation
and reproduction manifests identify its source; no external evolution dependency
has been incorporated. This decision covers the LL reference, not production NLL.

Primary capability pages inspected during this stage:

| Candidate | Verified information | Decision |
| --- | --- | --- |
| [eMELA](https://github.com/gstagnit/eMELA) | Documents space-like lepton PDFs, Fortran numerical evolution and a C++ analytical endpoint wrapper | Remains a candidate for NLL infrastructure; this information does not establish a ready time-like branch |
| [APFEL++](https://vbertone.github.io/apfelxx/html/namespaceapfel.html) | Documents time-like QCD initialisers and QED utilities | Remains a candidate; a complete time-like NLO QED matrix/interface has not been established |
| Standalone scalar LL | Existing LO kernel can be exponentiated and checked with analytic moments and independent real-space convolution | Selected for a small reference implementation without coupling production FF correctness to a PDF extension |

No repository has been forked or software revision pinned for eMELA/APFEL++.
The table is a documentation-level audit, not a claim of exhaustive source
inspection or a demonstrated build. A production NLL choice still requires a
source/build spike, a pinned revision and licence record if reused, verified
time-like kernel content, matrix orientation, endpoint treatment, threshold
support and fixed-order re-expansion. The standalone LL implementation does not
close the plan's NLL backend gate.

Version 0.8 adds a standalone coupled real-moment NLO reference, using the
complete one-family time-like QED matrix and two-loop MS̄ coupling. Its small
RK4 implementation is checked against independent continuation kernels,
analytic moments/running and Magnus matrix integration. This supplies a
validated reference for the next endpoint-controlled x-space implementation;
it does not select a production inversion/backend or implement thresholds.
See [moment evolution](moment_evolution.md). No external evolution or HPL
source has been incorporated.

Version 0.9 extends the standalone prototype to a signed energy-weighted
log-loss lattice, sharing the FFT and matrix integrator with the prior
references. Endpoint projection and weighted wrap/tail controls are tested
against real moments and an independent causal power-series solver. This is
the selected one-family x-space reference; production thresholds, observable
pair/hard matching and broader measurement validation remain separate gates.
See [x-space evolution](xspace_evolution.md).
