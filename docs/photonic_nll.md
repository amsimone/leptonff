# Fixed-alpha photonic NLL with unitary decay normalisation

The later [complete Z mass-radiator audit](z_radiator_validation.md) checks
the full hard-matched second coefficient independently. It also records
fixed-Q KKMC comparisons and the remaining physical validation limits.
The subsequent [McMule fixed-order benchmark](mcmule_fixed_order.md) adds
external O(alpha^2) moment evidence for the retained mass-log coefficients,
while resolving finite contributions beyond the stated NLL prescription.

Version 0.24 adds a separate `--photonic-nll` W/Z mode. It retains fixed
alpha, excludes pair channels and fermion-loop running, and includes the
time-like NLO photonic kernel. The existing `--photonic-ll` prediction is
unchanged. NLL within the photonic sector does not require adding pairs or
running; those define different physical-content choices.

## Source and the normalisation convention

Let `a=alpha/(2*pi)`, `L=log(mu_frag^2/mu0^2)` and
`U=exp(a*L*P0)`. We use the diagram-class photonic distribution `G` in
[Arbuzov–Melnikov, hep-ph/0205172, eq. (15)](https://arxiv.org/pdf/hep-ph/0205172)
and independently [Arbuzov–Scherbakova, 0706.2984, eq. (A.2)](https://arxiv.org/pdf/0706.2984).
The existing implementation and explicit delta coefficient are preserved.
As checked in the [kernel audit](nlo_kernel_audit.md),

```
g1 = G(N=1) = -13/4 + 3*zeta2 - 2*zeta3 = -0.719311605774...
```

This source-separated distribution is not a conserving scalar probability
kernel. Its first moment cancels against interference terms in the
pair-inclusive calculation; see hep-ph/0206036, eqs. (5.16)–(5.17). We do
not silently substitute a valence kernel containing pair interference, or
modify the source delta term to call it a different MSbar kernel.

The decay benchmark already fixes its inclusive normalisation through the
chosen branching ratio. Extend that explicitly to the photonic NLL shape:

```
E_raw(N) = U(N) * [1 + a*d1(N) + a^2*L*G(N)]
E_raw(1) = 1 + a^2*L*g1                 [to NLL]
E_unitary(N) = E_raw(N)/E_raw(1)        [expanded to NLL]
             = U(N) * [1 + a*d1(N) + a^2*L*(G(N)-g1)]
K_unitary = G - g1*delta(1-z).
```

Here `d1(1)=P0(1)=0`. This is an explicit conditioning of the source-sector
shape, not the unnormalised MSbar fragmentation function. The removed raw
factor and `g1` are exported in the callback log and run manifest. It does
not supply an absolute order-alpha-squared partial width. Conditional
normalisation is a theory convention, separate from the small numerical
number residual; no weights are clipped or fitted to PHOTOS.

## Matching and logarithmic counting

The implementation retains the entire LL tower and exactly one subleading
insertion globally: either `a*d1`, `a*H1`, or `a^2*L*K_unitary`.

```
W: U + a*U*d1 + a^2*L*U*K_unitary + a*H1[U]
Z: U x U
   + a*((U*d1) x U + U x (U*d1))
   + a^2*L*((U*K_unitary) x U + U x (U*K_unitary))
   + a*H1[U x U].
```

The original continuum order-alpha restoration and exact lattice overlap
subtraction remain in place. The existing finite W/Z currents and geometric
remainder are retained. NLO-kernel insertions on H1 or d1, two NLO-kernel
insertions, and two d1 insertions are beyond NLL and excluded. Thus hard
radiation still uses LL fractions; this is deliberate logarithmic counting.

For one leg the exported FF-only second coefficient is
`L^2*P0*P0/2 + L*(P0*d1 + K_unitary)`. The Z coefficient also includes the
two-leg first-order cross term with `d1 x d1` removed. In the NLL adapter
`ff_order2`/`sp_*_ff2` include this boundary/kernel coefficient. They are
explicitly FF-only, not the full hard-matched second coefficient. The LL
mode retains its historical auxiliary column definition.

## Numerics and validation

`PhotonicNLL` projects the signed resolved kernel onto log-loss cells and
transports it with the existing FFT LL operator. Unresolved contributions
are assigned to the atom according to the derived whole-plus conditioned
kernel. The signed convolution has explicit number and imaginary residuals.
A tilted total-variation bound controls periodic wrap. For log loss `v>=20`,
the source density satisfies the conservative bound
`|exp(-v)*G_regular(exp(-v))| <= 20*exp(-v)*(1+v^2)`; integrating it gives
the exported omitted-tail bound, including the factor two for replacing
the tail by an endpoint contribution. Domains failing either bound are
rejected. Cell and soft-cutoff errors still require refinement.

The dedicated C++ test checks the unchanged source, conditioned number,
literal signed convolution, one- and two-leg laboratory actions, unchanged
order-alpha overlap, second-order coefficients, cubic re-expansion,
LL-only H1 dressing, zero coupling and zero evolution. The independent
Python reference uses fixed transformed quadrature and an integral-defined
dilogarithm to check continuum moments and refinement through 65,536 cells.
`results/photonic_nll_reference_v024.json` records 55 checks.

The compiled adapter callback test uses four phase-space points for each
W charge and Z. It checks repeatability, all bin sums, unchanged Born/NLO,
hard and PHOTOS outputs, and active NLL shifts. These are implementation
checks, not an NLO-QCD integration or a precision certification.

The [v0.24 validation record](../benchmarks/photonic_nll_reference.md)
also records the frozen-LL callback regression and minimum-call W+/W-/Z
host smoke checks. One low-statistics W+ convergence warning is retained.

## Use

Build the current core and an isolated adapter with `build_adapter.py`.
The launcher rejects adapters without `photonic-nll-unitary-v1` capability
metadata and requires the applied NLL diagnostics in the run log.

```sh
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-photonic-nll-v024/mcfm_leptonff \
  --output results/new_photonic_nll_wplus --threads 1 \
  --photonic-nll --dress-hard --hard-samples 16 --spectra
```

The same option supports `wminus` and off-shell `z`. The comparison runner,
analyser and plotter recognise the distinct model `photonic_nll`:

```sh
python3 adapters/mcfm/run_comparison.py /path/to/MCFM-10.3.1 \
  --adapter build-mcfm-photonic-nll-v024/mcfm_leptonff \
  --results results/new_four_model_comparison \
  --models photonic photonic_nll leptons qed_hvp
python3 adapters/mcfm/analyse_comparison.py results/new_four_model_comparison \
  --models photonic photonic_nll leptons qed_hvp \
  --output results/new_four_model_comparison/report.json
```

These are reproduction commands, not runs performed for this increment.
The frozen v0.23 production results and plots do not contain this new curve.
The next production study should first address reusable host records and
configurable binning. Numerical/grid/cutoff and statistical precision remain
separate requirements. This mode is factorised decay at leading power in
the muon mass, not full NLO electroweak or mixed QCD–EW production.
