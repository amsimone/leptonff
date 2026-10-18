# Photon-parent hard matching and inclusive pair moments

Version 0.10 supplies `H1[gamma,mu±]`, its fiducial action at fixed fragmentation
fractions, and a complete hard-matched **pair-mass moment reference**. The
coupled fiducial x-space hard convolution remains pending. The model is a
pure unpolarised vector current, massless Born normalisation, one lepton family.
It is not a physical DY calculation.

## Pair counting

Sum inclusive identified opposite-charge pairs from distinct hard parents,
one observed lepton per selected direction. The Born parents are `(mu-,mu+)`.
The real three-body process supplies `(mu-,mu+)`, `(gamma,mu-)`, `(gamma,mu+)`.
Each photon-primary coefficient refers to **one specified primary charge**;
the photon FF may supply either observed charge. There is no extra-pair veto.

Two single-particle FFs do not supply correlations between two leptons from one
parent. A positive observed pair-mass cut excludes their unresolved collinear
limit. The moment weight `(m_pair²/M²)^(N-1)`, `N>=2`, also vanishes there.
Resolved same-parent contributions beyond the supplied hard order are outside
this reference's NLL coverage. The unweighted count at `N=1` is not supplied.

## Photon-primary coefficient

With `a=alpha/(2*pi)`, use massless invariants

```text
y = 2*k_gamma.k_selected_mu/M²,  t = 2*k_gamma.k_other_mu/M²,
q = 1-y-t,
x_gamma = y+t,  x_selected_mu = 1-t,  x_other_mu = 1-y,
R(y,t) = ((1-y)²+(1-t)²)/(y*t),  0<t<1-y.
```

R is the coefficient of a relative to massless vector Born. It follows by
changing variables in the existing real current, already checked against
Dirac traces. The massive scalar current is sourced to
[Bystritskiy et al., hep-ph/0505236v2, eq. (3)](https://arxiv.org/html/hep-ph/0505236v2).
At fixed y>0, only t=0 is collinear singular: the *unobserved* muon is parallel
to the photon, while the selected photon and muon remain separated.

For `r=m²/M²`, the massive t limits are

```text
t± = y*(1-y-2r ± sqrt((1-y)²-4r))/(2*(y+r)).
```

Integrating the exact massive current with its double-pole term gives

```text
g(y) = (1+(1-y)²)/y,
S_LP(y) = g(y)*(log(M²/m²)+2*log(1-y)-log(y))
          -4*(1-y)/y+(1-y)²/(2*y).
```

The double pole contributes `-2*(1-y)/y`; omitting it changes the finite term.
Subtract the Born `mu -> gamma` FF in MSbar:

```text
D1[mu -> gamma](y) = g(y)*(log(mu_frag²/m²)-2*log(y)-1),
H_gamma_mu(y) = g(y)*(log(M²/mu_frag²)+2*log(1-y)+log(y)+1)
               -4*(1-y)/y+(1-y)²/(2*y).
```

The FF is from
[Frixione, 1909.03886v2, section 6.3, eq. (293)](https://arxiv.org/html/1909.03886v2#S6.SS3).
The hard formula is our mass-regulator derivation, not a transcription from
that paper. It has no plus term or endpoint delta. The y=1 logarithm is
integrable; the y=0 number integral diverges. At mu_frag=M,
`integral dy y H_gamma_mu(y) = -13/3`.

`massive_vector_photon_pair_spectrum` integrates the exact current analytically
in t, using `y=(s_gamma_mu-m²)/M²`, with endpoint `1-2m/M`. Tests integrate the
independent Dirac current at fixed y and verify decreasing mass corrections.
The sum of both charge-specific y moments also reproduces the photon-energy
moment of the existing massive pair-mass spectrum. This is internal independent
current validation, not an external-generator comparison.

## Fixed-fraction fiducial action

At fragmentation fractions `(fg,fm)`, the observed pair mass is
`M*sqrt(y*fg*fm)`. Energies are `M*(y+t)*fg/2` and `M*(1-t)*fm/2`, with opening
cosine `1-2*y/((y+t)*(1-t))`. Common-rotation band acceptance gives A(y,t).
The collinear projection is `B(y)=isotropic_acceptance(M,y*fg,fm,cuts)`, so

```text
H_gamma_mu[F] = integral dy H_gamma_mu(y)*B(y)
             + integral dy dt R(y,t)*(A(y,t)-B(y)).
```

The second integral is locally finite. `fiducial_vector_photon_pair_hard`
requires mass_min>0, integrates only actual mass support, and exposes quadrature
and t/(1-y) sliver controls. `HardDecayOptions.soft_cutoff` is unused: the positive
mass lower cut excludes soft photons exactly. Sharp angular cuts can slow sliver
convergence; the reference refinement is not a universal error bound.

The API output includes neither a nor photon-conversion FFs nor the sum over
primary charges. It is not itself a muon-pair rate.

## Matched inclusive moments and coupling

`vector_pair_moment` joins all hard-parent channels. Write
`S=D[mu- -> mu-](N)`, `O=D[mu- -> mu+](N)`, and
`G=D[gamma -> mu-](N)=D[gamma -> mu+](N)`. With Hm the existing muon-pair hard
mass moment and Hg the photon-primary hard moment,

```text
R_N = (S²+O²)*(1+a_h*Hm) + 2*a_h*Hg*G*(S+O),
a_h = alpha_MSbar(M)/(2*pi).
```

The factor two sums primary charges; both observed charge assignments are
already included in G*(S+O). The photon term starts at a²*log(M²/m²), so it
is required at NLL. These N>=2 moments are angular inclusive; the driver does
not apply the fiducial cuts in the configuration.

One explicit MSbar input at mu_0 determines FF running and a_h. The hard
coupling stays at M when varying mu_frag. The default 1/137.035999084 is a model
MSbar input, not an on-shell conversion. Two-loop running rows supply the
stated one-family NLL reference; fixed/one-loop rows are running diagnostics.
No physical thresholds or hadronic vacuum polarisation are supplied. Products
and exponents retain some terms beyond NLL. Constant two-loop FF boundaries
and H2 are omitted.

The API also exposes the expansion in a0=alpha_MSbar(mu_0)/(2*pi). For
`D=I+a0*C1+a0²*C2+...` and `Lh=log(M²/mu_0²)`,

```text
R1 = 2*C1_same + Hm,
R2 = 2*C2_same+C1_same²+C1_opposite²
     +Hm*(2*C1_same+beta0*Lh)+2*Hg*C1_gamma_to_mu.
```

beta0=2/3 for running modes, zero for the fixed diagnostic. Tests recover the
complete O(alpha) leading-power massive pair coefficient and a cubic remainder
after R1/R2 subtraction. Varying mu_frag leaves R1 invariant; the residual R2
scale shift has no mass logarithm. Omitting the photon channel or using a0
for the hard coupling spoils that NLL check. Nonzero full-result scale shifts
include unsupplied NNLL terms and are not total physics-error estimates.

Python reconstructs FFs with analytic running/Magnus steps, contracts unreduced
parent/observed-charge tensors and checks hard moments with harmonic sums.
The separate fiducial check integrates allowed azimuth intervals, and recovers
the projected coefficient by subtracting the FF from a numerically integrated
massive current with a shrinking auxiliary mass. See the
[benchmark](../benchmarks/pair_matching_reference.md) for achieved errors.
