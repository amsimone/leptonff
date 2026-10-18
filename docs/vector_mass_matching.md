# Angular-inclusive vector-current hard matching

Version 0.4 adds a restricted first-order hard-decay benchmark. The parent is
at rest and unpolarised, its coupling is purely vector, and the only cut is
`mass_min < m_mumu < mass_max`. Rates are relative to **massless vector Born**,
not a fitted Z width. The exact massive Born ratio is beta*(1+2*r), with
r=m_mu^2/M^2 and beta=sqrt(1-4*r). This module does not apply pT/eta cuts,
include axial mass effects, or supply the fully differential hard Z decay.

## Factorisation and finite coefficient

Write a=alpha/(2*pi), z=m_mumu^2/M^2, p(z)=(1+z^2)/(1-z), and use whole-kernel
plus prescriptions throughout. Expanding the massive real radiator at fixed
0<z<1 and including the inclusive correction gives

```text
C1_LP = (3/2) delta(1-z) + 2 [p(z)*(log(M^2/m_mu^2)+log(z)-1)]+ .
2 D1  = 2 [p(z)*(log(mu_frag^2/m_mu^2)-2*log(1-z)-1)]+ .
H1    = (3/2) delta(1-z)
        + 2 [p(z)*(log(M^2/mu_frag^2)+log(z)+2*log(1-z))]+ .
```

The two FF insertions give 2*D1 since either leg alone changes the pair mass
at this order. H1 is finite and mass independent. The implementation calculates
H1+2*D1 separately and checks it against C1_LP, including cancellation of
mu_frag and mu_0 dependence. The inclusive coefficient is 3/2; the N=2 hard
moment at mu_frag=M is 61/6. No hard term is inserted into the existing LL or
fiducial drivers by this module.

## Massive reference and source audit

The real spectrum uses [Bystritskiy et al., hep-ph/0505236v2](https://arxiv.org/html/hep-ph/0505236v2),
eq. (5), with photon fraction nu=1-z. Its prefactor divided by massless Born
is a. Cumulatives use numerical integration with a quadratic threshold map;
real-only bins are integrated directly. The inclusive coefficient uses
[Czyz et al., hep-ph/0404078v2](https://arxiv.org/html/hep-ph/0404078v2), eq. (19),
converted from alpha/pi and massive Born to this module's conventions.

The latter source's eq. (17), as displayed in arXiv HTML and PDF text extraction,
gave a discrepancy in its `(1+beta_m^2)` logarithmic term. Changing that factor
to `(1+beta_mu^2)` agreed with the independent current, but this is a local
inference, not a confirmed erratum. For example, at M=91.1876 GeV, m=1 GeV,
z=0.7, the extracted formula's derivative differed from the trace by 0.07655
in dC/dz. The implementation avoids eq. (17) entirely; it uses the independent
differential formula above. No author confirmation is claimed.

## Independent checks

`vector_current.cpp` constructs the two-leg photon-emission current using
explicit Dirac matrices and on-shell spin projectors. It sums two transverse
photon polarisations and three parent spin states. The factor multiplying the
trace is `(1-z)*beta_pair/16`, after division by massless vector Born. A
hyperbolic angular map resolves the mass-regulated collinear peaks. Tests check
on-shell kinematics, the photon Ward identity, and the angular integral against
the scalar spectrum for four masses and three z points, including a large mass.

Python uses Simpson integration, a separate logarithmic soft map, the integral
definition of Li2, and virtual-plus-soft eq. (13) of hep-ph/0505236v2. Combining
that term with real radiation at two soft separators and extrapolating linearly
tests the C++ inclusive coefficient without reusing its formula. The manifest
records the separator shifts; their extrapolated residual is not a Monte Carlo
uncertainty. Separate inclusive, real-tail and mass-scan fixtures check both
branches of the mass-window action. The upper mass cut exactly at M is rejected:
with strict cuts it would exclude Born while accepting arbitrarily soft radiation.

This established a mass-only part of G2 in version 0.4. Version 0.5 adds
[symmetric fiducial cuts and additive NLO+LL](fiducial_vector_matching.md).
Version 0.6 adds [LL dressing of H1](hard_dressing.md). The NLL kernel,
physical spin-correlated hard subtraction and Z-current validation remain
required for the final matched-decay gate.
