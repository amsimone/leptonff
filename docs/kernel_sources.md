# Kernel provenance

Source: [Frixione, arXiv:1909.03886v2, section 6](https://arxiv.org/html/1909.03886v2#S6).
Set finite scheme functions to zero (MS̄), unit charge, one muon family, and
`L0 = log(mu_0²/m_mu²)`. Code stores coefficients of `alpha/(2*pi)`.

| Parent → observed | O(α) coefficient | Source equations |
| --- | --- | --- |
| μ± → μ± | `[(1+z²)/(1-z) * (L0 - 2 log(1-z) - 1)]+` | (289), (142) |
| μ± → γ | `(1+(1-z)²)/z * (L0 - 2 log(z) - 1)` | (293), (210) |
| γ → μ± | `(z²+(1-z)²) * L0` per charge | (300), (255) |
| γ → γ | `-(2/3) L0 delta(1-z)` | (302), (260) |
| μ± → μ∓ | zero at this order | order counting |

LO kernels are their derivatives with respect to `L0`. Momentum sums include
both photon-conversion charges. The separate NLO time-like charged-lepton
block is implemented in `nlo_lepton.hpp`. Version 0.8 adds the full one-family
matrix through `nlo_qed.hpp` and running real-moment evolution. The
[photon derivation](nlo_photon_kernels.md) and [solver contract](moment_evolution.md)
record normalisations and checks. Threshold and fiducial NLL matching remain open.

Independent analytic diagonal moment checks used in tests:

```text
d1(N=1) = 0
d1(N=2) = -(4/3)L0 - 17/9
d1(N=3) = -(25/12)L0 - 265/72
```

The non-singlet LL solver exponentiates this same LO diagonal whole-plus kernel.
See [LL evolution](ll_evolution.md) for the Mellin solution, scalar-theory scope
and regulator treatment; it does not import NLO space-like kernels.
The signed diagonal boundary is now transported by that LL operator; see
[boundary matching](boundary_matching.md). The incomplete
[NLO source audit](nlo_kernel_audit.md) records the direct QED lepton-block
transcription, pair-interference charge check and complete one-family matrix.
