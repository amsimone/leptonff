# Coupled NLO time-like moment evolution

Version 0.8 supplies a standalone real-Mellin-moment reference for the complete
one-family MS̄ matrix. `evolve_qed_moment(N,parameters,options)` supports real
`N>1` and returns the propagator, boundary-transported FFs, running-coupling
integrals and explicit first/second expansion coefficients. Matrix rows are
parents, columns are observed species, in order `{mu-,mu+,gamma}`.

## Theory and boundary contract

The input is named `alpha_msbar_at_mu0` explicitly. Its default numerical
value is a model input, not a computed conversion of the on-shell coupling.
There is exactly one active unit-charge lepton family at every scale and no
threshold or hadronic contribution. The [NLO matrix](nlo_photon_kernels.md)
includes inclusive pairs and both charges. Fixed and one-loop running modes
are diagnostics; two-loop mode supplies the running used for the NLL reference.

With `a=alpha_MSbar/(2*pi)`, `L=log(mu_frag²/mu_0²)` and `E(0)=I`, solve

```text
da/dL = b0 a² + b1 a³,       b0=2/3, b1=1  (two-loop mode)
dE/dL = (a P0 + a² P1) E
D(N)  = E(N) (I + a_initial d1(N,mu_0,m_mu))
```

The beta coefficients follow from
[Bertone et al., arXiv:1911.12040v3](https://arxiv.org/pdf/1911.12040),
equations (4.11)–(4.12), after converting from alpha to `a`. This uses their
coupling convention, not their space-like kernels. The FF boundaries are the
existing MS̄ one-family coefficients in [kernel sources](kernel_sources.md).
The boundary is applied once, on the right of the evolution matrix.

The NLO-generator solution contains the NLL towers and some terms beyond NLL;
it is not a globally linear insertion prescription. Unknown constant two-loop
boundaries are omitted. This is a separate reference from the existing LL
fiducial driver, which retains at most one finite hard/boundary insertion.
The two prescriptions must not be combined by substituting a new alpha value.

The full matrix does not have finite N=1 entries because of singlet/survival
small-z singularities. The API rejects N<=1; the finite charge-odd N=1 kernel
is checked separately. At admissible N the charge-odd sector agrees with
`exp(P0_valence integral(a) + P1_valence integral(a²))`.

## Expansion and integration checks

The returned coefficients are in powers of the **initial** `a`:

```text
D = I + a_initial (d1 + L P0)
      + a_initial² [L P1 + L P0 d1 + (L²/2)(P0 P0 + b0 P0)] + ...
```

For fixed coupling set `b0=0`. The `L P1` term is explicit and distinct from
the earlier boundary-generated `L P0 d1`. The two-loop beta coefficient first
affects higher orders; it is included in the resummed solution. A coupling
scan verifies cubic onset of the residual after subtracting both coefficients.

C++ integrates the coupling, its two moments and all matrix entries jointly
with RK4, doubling the step count until the maximum absolute state shift meets
`step_tolerance`. Kernel quadrature has its own integration controls. Neither
step refinement nor moment sum rules bound an x-space inversion error.
No entry is clipped or renormalised to force a conservation identity.

Running matrices at different scales need not commute. Tests verify
`E(mu2,mu0)=E(mu2,mu1) E(mu1,mu0)` using the transported MS̄ coupling,
including a larger-coupling fixture that resolves matrix-ordering errors.
Backward propagation gives the inverse. Equal scales give the identity
propagator and retain the boundary; zero coupling gives the identity FF.
Unsupported strong running and failed refinement throw.

The independent Python check obtains the coupling from its analytic implicit
two-loop equation. It integrates the matrix in `s=log(a/a_initial)` using
fourth-order Magnus steps and matrix exponentials, including the commutator.
LO/boundary moments are analytic and NLO photon moments come from the separate
space-like-plus-continuation derivation. It checks the full matrices and both
expansion coefficients, not just a diagonal observable.

## Remaining observable work

This reference supplies real moments, not complex-N inversion or an x-space
FF action. Version 0.9 supplies an independently checked
[energy-weighted coupled grid](xspace_evolution.md), including endpoint and
domain controls, moment recovery and a Born FF fiducial reference. Before full
fiducial NLL use, align
the renormalised coupling with the hard calculation, define inclusive versus
exclusive extra-pair selection and count photon hard-parent channels.
Physical flavour thresholds and spin-correlated boosted decay remain open.
See [results and reproduction](../benchmarks/moment_reference.md).
