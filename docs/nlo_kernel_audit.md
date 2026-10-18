# Time-like NLO kernel audit — lepton block implemented

Version 0.24 uses the source photonic class in a separately derived
[unitary-conditioned photonic NLL shape](photonic_nll.md). Its inclusive
source factor is retained as metadata; the source distribution below is
unchanged. This is not a substitution of G for the conserving valence kernel.

Version 0.7 adds the one-family, unit-charge MS̄ NLO charged-lepton block,
including pair terms, in `nlo_lepton.hpp`. The coefficient multiplies
`a²=(alpha/(2*pi))²` in evolution. The LL decay drivers do not use it yet.
Version 0.8 completes the [photon channels](nlo_photon_kernels.md) and adds
[coupled real-moment evolution](moment_evolution.md) with two-loop running.
Version 0.9 adds the [coupled x-space reference](xspace_evolution.md).
Threshold and full fiducial hard/pair matching remain open.

## Direct QED transcription

[Arbuzov–Melnikov, hep-ph/0205172v1](https://arxiv.org/pdf/hep-ph/0205172),
equations (14)–(19), supplies the same-charge photonic, non-singlet-pair,
singlet-pair and interference terms. The photonic formula also appears in
[Arbuzov–Scherbakova, arXiv:0706.2984v1](https://arxiv.org/pdf/0706.2984),
appendix (A.2). Write `p=(1+z²)/(1-z)`, `l=log(z)`, `s=log(1-z)` and
`c=3/8-3*zeta2+6*zeta3`. The implementation transcribes:

```text
G = c delta + p (2 l s - 2 l² - 2 Li2(1-z))
              + (1+z) l²/2 + 2 z l - 3 z + 2
B = (-4 zeta2/3 - 1/6) delta - (20/9) [1/(1-z)]+
              - (2/3) p l - 2/9 + 22 z/9
S = (1+z) l² + (-5 - 9 z - 8 z²/3) l
              - 8 - 20/(9 z) + 4 z + 56 z²/9
I = p (2 Li2(1-z) + 3 l/2) - (7/2)(1+z) l - 7 + 8 z
```

Only the displayed `1/(1-z)` term carries a plus prescription. No endpoint
renormalisation is added. The code uses the standard positive `Li2(x)` for
`0<x<1`, explicitly defined in the 2002 paper's equation (19).

[Arbuzov, hep-ph/0206036v2](https://arxiv.org/pdf/hep-ph/0206036), equation
(5.15), supplies the opposite-charge interference contribution:

```text
J = 2 (1+z²)/(1+z) (-2 Li2(-z) - 2 l log(1+z) + l²/2 - zeta2)
              + 2 (1+z) l + 4 (1-z)
```

The reconstructed parent→observed lepton block is
`P1[same]=G+B+S+I`, `P1[opposite]=S+J`; charge conjugation fixes both
parent signs. These are inclusive identified-lepton densities with extra
pairs, not an exclusive two-muon selection.

For an independent density comparison, adding `G+I` cancels `Li2(1-z)`.
The resulting `G+B+I` agrees at `0<z<1` with the non-singlet time-like
expression in [Blümlein–Kawamura, hep-ph/0211191v1](https://arxiv.org/pdf/hep-ph/0211191),
equation (30). This checks the regular density; endpoint coefficients come
explicitly from the first source.

## Conservation and executable checks

The 2003 paper's equations (5.16)–(5.21) explain the first-moment bookkeeping.
Let `k=13/4-3*zeta2+2*zeta3=0.719311605774509...`. Our integrals give:

```text
G(N=1) = -k       B(N=1) = 0
I(N=1) = 2 k      J(N=1) = k
V = P1[same] - P1[opposite] = G+B+I-J
V(N=1) = 0
```

Thus the nonzero photonic first moment is **not** repaired by adding a
whole-plus prescription or changing its delta coefficient. Also, `I=2*J`
holds for the first moments only, not as functions. The source-separated
photonic term must not be substituted for a conserving scalar kernel merely
because the existing LL backend is called photonic/non-singlet.

`S` has a `1/z` singularity. Both inclusive charge channels require `N>1`;
requesting their number moments throws. `nlo_valence_kernel()` cancels `S`
analytically before distribution construction and has a finite `N=1` moment.
Subtracting two already-divergent number integrals is never used.

`leptonff_kernels` reports all five components, the assembled lepton block,
valence moments at `N=1,2,3,5` where finite, a finite-bound action and a
tolerance refinement. The Python checker uses integral-defined dilogarithms
and endpoint-transformed fixed Gauss quadrature independently of C++'s
series/adaptive integrator. Analytic harmonic-sum moments check `B`, `S` and
the log-only `G+B+I` form. See the
[numerical benchmark](../benchmarks/nlo_kernel_reference.md).

The original `nlo_lepton_kernel` API remains limited to charged leptons and
throws for photon parents or daughters. The new `next_to_leading_kernel`
supplies the complete matrix, including photon survival; every parent momentum
sum passes. The coupled real-moment reference includes two-loop running and
O(α²) extraction. The x-space reference now recovers those moments and supplies
FF-only fiducial matching; full hard/pair matching remains required before a
fiducial NLL result.

## Earlier QCD/backend audit

[Stratmann–Vogelsang, hep-ph/9612250v1](https://arxiv.org/html/hep-ph/9612250)
section 3 covers unpolarized NLO time-like evolution. Equation (18) uses
`alpha_s/(2*pi)`, matching this project's expansion convention. Equations
(34)–(36) give analytic-continuation relations and their corrections, **not
standalone complete NLO kernels**. The complete unpolarized formulas are referred
back to the earlier calculations. The text after (39) specifies endpoint
information and an existing sign misprint. These details must accompany a
transcription; the later polarized section is not the target kernel set.

[Mitov–Moch–Vogt, hep-ph/0604053](https://arxiv.org/html/hep-ph/0604053)
uses `alpha_s/(4*pi)` and provides higher-order non-singlet continuation checks.
For coefficient tables in that convention, conversion to the project's `a`
requires `P0_here=P0_there/2` and `P1_here=P1_there/4`, before any charge/colour
mapping. It is a cross-check, not a complete QED matrix source.

The [APFEL++ time-like header documentation](https://vbertone.github.io/apfelxx/html/splittingfunctionsunp__tl_8h.html)
identifies QCD classes `P1Tnsp`, `P1Tnsm`, `P1Tps`, `P1Tqg`, `P1Tgq`, `P1Tgg`.
They are concrete candidates for a source/build comparison. Their existence
does not establish a complete time-like QED backend or fix this project's
parent/daughter index orientation. The repository's [licence file](https://github.com/vbertone/apfelxx/blob/master/LICENSE)
contains GPL version 3. No APFEL++ code has been incorporated into this project.

## Required channel reconstruction

The lepton rows have a direct QED transcription above. Version 0.8 reconstructs
the photon rows independently through full analytic continuation with explicit
per-charge and closed-loop normalisations.

| Project channels | QCD information needed | Checks before acceptance |
| --- | --- | --- |
| μ± → μ± | Implemented direct QED `G+B+S+I` | Explicit delta/plus terms; independent moments and density comparison |
| μ± → μ∓ | Implemented direct QED `S+J` | Per-charge singlet; interference is not a flavour-count substitution |
| μ± → γ | Implemented time-like radiation channel | LO mapping, independent continuation, momentum and finite-bound checks |
| γ → μ± | Implemented conversion, per charge | Separate produced charge sum and closed-loop multiplicity; both signs checked |
| γ → γ | Implemented real and virtual contributions | Explicit delta=-1; per-parent momentum; one-family beta content |
| Other lepton families | Pure singlet/flavour transfer and their thresholds | Explicit flavour basis and charge factors; no inference from the one-family diagonal alone |

The original unpolarized sources identified in these references are Curci,
Furmanski–Petronzio, Nucl. Phys. B175 (1980) 27; Furmanski–Petronzio,
Phys. Lett. 97B (1980) 437; and Floratos–Kounnas–Lacaze, Nucl. Phys. B192
(1981) 417. A formula-level transcription from these or an independently
validated modern implementation can provide further cross-checks. The photon
derivation and authors' reference-routine audit are now recorded separately.

## Acceptance procedure

1. Pin the source version and record every regular, plus and delta coefficient
   separately, including errata and the expansion normalisation.
2. Fix parent/daughter and NS±/singlet conventions by reproducing the existing
   LO matrix first. State the exact real/virtual fermion content.
3. Abelianise by charge-labelled graph classes. Do not assume that selecting
   one colour monomial in an inclusive QCD result removes every pair term.
4. Check charge conservation in the reconstructed valence channel and the
   complete per-parent N=2 momentum sums, including photon survival.
5. Compare selected Mellin moments and x-space values with an independent
   calculation/backend and verify the time-like/space-like difference.
6. Re-expand the evolution to O(α²), checking the **additional** `L P1` term
   separately from the already implemented `L P0⊗d1` and `L² P0⊗P0/2` terms.

Version 0.8 passes these checks for real moments in the one-family model.
Version 0.9 adds an endpoint-controlled x-space action and independent causal
checks. Production NLL still requires consistent hard coupling, all retained
hard-parent channels and extra-pair matching; FF-only tests do not close those
observable gates.
