# One-family NLO time-like photon channels

Version 0.8 completes the unit-charge MS̄ `{mu-,mu+,gamma}` NLO matrix in
`next_to_leading_kernel(parent,observed)`. It uses `a=alpha/(2*pi)` and
`mu_ren=mu_frag`. The lepton block retains the direct QED transcription in
the [earlier audit](nlo_kernel_audit.md). The photon expressions below are our
derivation and numerical cross-check, not a claim that a PDF kernel is already
a fragmentation kernel.

## Independent source routes

The starting space-like QED expressions are
[de Florian–Sborlini–Rodrigo, arXiv:1606.02887v3](https://arxiv.org/pdf/1606.02887v3),
equations (55), (61), (62) and (68). Set one unit-charge lepton family and no
quarks; equation (61) removes the quark colour multiplicity.
Apply the **full** continuation, including its finite correction, from
[Stratmann–Vogelsang, hep-ph/9612250v1](https://arxiv.org/html/hep-ph/9612250),
equations (34), (36), (38), (39). In particular, crossing a squared endpoint
logarithm contributes `log²(1-x) -> log²(x-1)-pi²`.

For interior `z`, writing `l=log(z)`, `s=log(1-z)`,
`g=(1+(1-z)²)/z` and `p=z²+(1-z)²`, the per-charge continuation corrections
used by the independent checker are:

```text
C[mu->gamma] = 2 g (4 Li2(1-z) + 2 l s + s)
                + (-4 z (3-z) l - 20 + 28 z - 5 z²)/(2 z)
C[gamma->mu] = (4/3) p (1-2 l) + 2 p (-4 Li2(1-z) - 2 l s + s)
                + (-4 (2-z)(1-2z) l - 5 + 28 z - 20 z²)/2
C[gamma->gamma] = 4/(9 z) ((1-z)(38+47z+38z²)
                + 3(1+z)(4+11z+4z²) l)
```

The time-like regular parts are `z*AC(S[mu<-gamma](1/z))+C[mu->gamma]`,
`z*AC(S[gamma<-mu](1/z))+C[gamma->mu]` and
`-z*AC(S[gamma<-gamma](1/z))+C[gamma->gamma]`, respectively. Here the arrows
on `S` denote daughter←parent PDF indices. The Python implementation retains
these crossed expressions; C++ uses the algebraically simplified results.

As an additional audit, the coefficient functions in the authors'
[NLO time-like reference routines](https://www.liverpool.ac.uk/~avogt/xpij1te.f)
give the same results. They use `alpha_s/(4*pi)` and an implicit `TR=1/2`.
Restoring external charge multiplicities and closed-loop factors gives these
weights on their indicated colour polynomials **before dividing by four**:

| Physical channel | Source polynomial | One-family QED weight |
| --- | --- | ---: |
| mu→gamma | `X1GQTA: CF²` | 1 |
| gamma→each mu charge | `X1QGTA: nf*CF` | 1 |
| gamma→each mu charge, closed-loop term | `X1QGTA: nf²` | 2 |
| gamma→gamma regular and local | `X1GGTA/X1GGTC: nf*CF` | 2 |

The produced-charge sum and a closed Dirac loop are different factors. For
example, `2*nf*TR` for the observed charge sum becomes two and is divided by
two for each observed muon charge; a closed `TR*nf` becomes one. Terms with
non-Abelian vertices vanish. This channel-specific bookkeeping also reproduces
the existing LO photon conversion and radiation normalisations. A single global
`nf` substitution would miss the relative weights. No external Fortran or HPL
implementation has been incorporated; the fixed-version paper continuation is
the reproducible derivation, independent of changes to the authors' website.

## Implemented expressions and checks

The compact coefficients, all in this project's `a²` convention, are:

```text
P1[mu->gamma] = -1/2 + 9z/2 + (-8+z/2) l + 2 z s + (1-z/2) l²
                 + g (s² - 4 l s - 8 Li2(z))

P1[gamma->mu] = -6 + 23z/2 - 10z² + (-5/2+2z+2z²) l
                 + (-3+2z-2z²) s + (-1/2+z-2z²) l²
                 + p (-2 zeta2 + 8 Li2(z) + 6 l s - s²)
                 - 4/3 + p (-8/9 - (4/3)(l+s))

P1[gamma->gamma] = -delta(1-z) - 4 + 12z - 164z²/9 + 92/(9z)
                 + (10+14z+16z²/3+16/(3z)) l + 2(1+z) l²
```

There are no plus terms in these photon channels. The virtual endpoint is
`-delta(1-z)`, derived from the source, not chosen to enforce momentum after
integration. The off-diagonal NLO coefficients have finite `N=1` moments,
despite the divergent photon-number moment of the LO radiation kernel.
The NLO photon-survival coefficient requires `N>1`.

The complete matrix satisfies each parent momentum sum at N=2. The photon
entries give `P1[gamma->each mu](2)=65/54` and
`P1[gamma->gamma](2)=-65/27`. The muon radiation entry cancels the sum of
the already checked same/opposite-charge lepton entries. Moments at N=3,5,
finite-bound actions and independent endpoint refinements also pass.
The [benchmark](../benchmarks/moment_reference.md) records achieved residuals.
