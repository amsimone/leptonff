# Investigation of the fixed-Q KKMC residuals

The subsequent [strict fixed-order McMule study](mcmule_fixed_order.md) has now
performed the coefficient test proposed below. Its six-Q moments support the
implemented mass-logarithm content, with finite second-order remainders; the
full resummed residual remains incompletely classified. The results below
describe the preceding resummed coupling/mass investigation.

The residuals are reproducible, have a substantial mass dependence, and cannot
be accounted for by the tested finite second-order completion alone. Their
division into a second-order logarithmic mismatch and higher-order differences
is **not yet established**. There is no justified change to the production
kernel from these tests. In particular, neither a successful coefficient audit
of our prescription nor a statistically resolved generator difference establishes
the prescription's physical NLL accuracy.

## Controlled calculations

All new calculations use the original high-statistics KKMC executable/defaults
snapshots. The upstream source tree, commit and Fortran flags agree with the
earlier small-alpha controls; their driver update changed metadata and the
radiation-off branch. Every binary and raw-input hash is retained in the report.
All streams are independent. No changed coupling, mass, cutoff, or prescription
is pooled with a different physical point.

At Q=91.1876 GeV the additional controls are:

| Mass | alpha / alpha_physical | New weighted trials |
| --- | ---: | ---: |
| m_mu | 1/2 | 200 million |
| m_mu | 2 | 200 million |
| m_mu/100 | 1 | 200 million |
| m_mu/100 | 1/2 | 200 million |
| m_mu/100 | 2 | 50 million |

These add 850 million weighted trials. The analysis also uses the existing
200-million-trial physical point and the 20-million-trial alpha/8 and alpha/16
controls. This is Born QCD and fixed-alpha photonic FSR, with ISR, IFI, weak,
pairs and running disabled. These are theory diagnostics, not physical changes
of the muon mass or coupling for phenomenology.

Below, delta_N = R_N(CEEX2) - R_N(leptonff NLL), and ppm denotes 1e-6 in
the normalized moment itself, **not** a fractional change divided by R_N.
Errors are MC standard errors, with normalization/measurement covariance.

## Coupling dependence

For the N=2 moment at the physical muon mass:

| alpha factor f | delta_2 [ppm] |
| ---: | ---: |
| 1/2 | -20.30 +/- 4.90 |
| 1 | -86.06 +/- 7.13 |
| 2 | -441.11 +/- 10.71 |

The ratios to f^2 are -81.18, -86.06, -110.28 ppm. The residual is therefore
compatible with a leading quadratic effect plus higher orders; a purely
quadratic description of all five coupling points has chi-square 19.30 for
four degrees of freedom (p=0.00069). A fit delta=c2*f^2+c3*f^3 gives

    c2 = -63.24 +/- 12.71 ppm
    c3 = -23.51 +/-  6.72 ppm
    chi-square = 7.08 / 3 degrees of freedom.

These coefficients are diagnostics of that polynomial model. They are not
determinations of physical fixed-order coefficients. The older alpha/8 point
has a 2.63-sigma fit residual, correlated across the moments. Adding a linear
term gives c1=18.59 +/- 26.66 ppm, compatible with zero. Allowing a quartic
term makes the extrapolation much less stable: c1=154.08 +/- 72.74 ppm and
c2=-617.37 +/- 265.06 ppm. The data do not justify declaring a first-order
problem, or quoting the quadratic/cubic decomposition without its truncation
qualification. The report contains all fits and their absolute MC covariance;
fit covariance is not rescaled to force chi-square per degree of freedom to one.

## Mass dependence

At the physical coupling the N=2 residual changes from -86.06 +/- 7.13 ppm
at m_mu to -326.38 +/- 9.98 ppm at m_mu/100. Define a=alpha/(2*pi),
L=log(Q^2/m^2), and strip the common LL evolution U_N=exp(2*a*L*P_N):

    D_N = delta_N/(a^2*U_N).

For N=2, D changes from -66.53 +/- 5.51 to -259.61 +/- 7.94 as L changes
from 13.5209 to 22.7313. That is a 20-sigma change at fixed physical coupling.
The single finite-completion ansatz a^2*U_N*constant cannot describe it.
Finite-coupling mass dependence can arise from retained or omitted higher-order
logs as well as from a second-order logarithmic mismatch.

The extra mass dependence, delta(m/100)-delta(m), is:

| alpha factor f | Difference [ppm] |
| ---: | ---: |
| 1/2 | -41.68 +/- 8.38 |
| 1 | -240.32 +/- 12.26 |
| 2 | -1125.19 +/- 32.33 |

Fitting this difference to f^2 and f^3 gives a quadratic coefficient
-181.97 +/- 22.55 ppm. However, f^3 and f^4 also fit it acceptably. When
all three powers are allowed, the quadratic coefficient becomes
-57.74 +/- 92.76 ppm, with no remaining fit degrees of freedom. Thus this
scan does **not** uniquely identify a missing alpha^2*L term. There is only
one independent mass lever arm, so it cannot separate L from L^2 either.

## Source and convention checks

The library's regular time-like kernel and delta coefficient agree with
[Arbuzov--Melnikov, eq. (15)](https://arxiv.org/html/hep-ph/0205172) and
[Arbuzov--Scherbakova, appendix A.2](https://arxiv.org/html/0706.2984).
The latter distinguishes the space-like and time-like functions explicitly.
The source's photonic first moment is nonzero. Its conditioning must be
kept distinct from the pair-interference cancellation discussed in
[Arbuzov, section 5.2](https://arxiv.org/html/hep-ph/0206036v2).
The comparison provides no basis for substituting a space-like, valence,
or pair-inclusive kernel to improve numerical agreement.

The pinned KKMC source has these relevant features:

* CEEX constructs amplitudes with zero-, one- and two-hard-photon residuals,
  multiplied by the YFS soft factors. Soft exponentiation does not by itself
  identify its arbitrary-order collinear series with DGLAP LL/NLL evolution.
  Higher-order completion differences therefore remain a candidate. See
  [the CEEX formulation](https://arxiv.org/abs/hep-ph/0006359).
* The FSR one-loop correction is in `BVR_MakeVfin`, used by `GPS_HfinPlus`.
  Comparing it to [Jadach et al.](https://arxiv.org/html/hep-ph/0109279)
  requires consistent parent/pair invariants, YFS subtraction and mass terms.
  A different-looking log coefficient is not sufficient evidence of a bug.
  No isolated incorrect term has been demonstrated in that implementation.
* The paired CEEX2-minus-EEX2 difference has a resolved component linear in
  alpha. Its fitted N=2 coefficient at physical coupling is -388.11 +/- 5.14
  ppm. These two weights cannot be substituted for one another as if they
  differed only by the same fragmentation LL-to-NLL increment.
* Weight 74 is not a complete third-order FSR upgrade. In `QED3.f`, the new
  third-order terms enhance ISR; the FSR virtual factor stops at delf2, and
  there is no pure three-FSR-photon residual in its total. Consequently it
  cannot serve as a complete alpha^3 FSR control.

The inspected files are the pinned upstream `KK2f/CEEX.f`, `GPS.f`, `BVR.f`,
`QED3.f`, `KK2f.f` and `.KK2f_defaults`. Their source-tree identity is in
the diagnostic report's backend manifests.

There is also a technical support restriction: the default parameter 638,
MinMassCEEX, is 1 GeV. `KK2f_Make` leaves CEEX weights zero below this pair
mass. Its normalized denominator therefore excludes that region, whereas the
FF reference covers its leading-power full support. Integrating the known
leading-power first-order real spectrum over the omitted physical region
estimates missing fractions of 0.70 ppm at m_mu and 3.27 ppm at m_mu/100.
Removing a low-z tail raises the normalized N>1 moments, opposite to the
observed deficit. This first-order effect is too small to explain the measured
residual, but the estimate is **not** an all-order bound or an exact threshold
calculation. An identical support definition is still required for a definitive
comparison. The driver now exports `ceex_min_mass`; pooling rejects different
support floors, and the comparison warning states the legacy restriction.
A separate 10,000-trial metadata smoke test passes.
Using the same seed in the original and updated drivers gives bitwise identical
observable means and covariance matrices for all 10,000 trials.

## The decisive next calculation

A fixed-order alpha^2 shape benchmark is more decisive than another fit to
independent resummed MC samples. It is also simpler than computing an absolute
two-loop width. Write the inclusive and measured rates as

    sigma = sigma0 + a*sigma1 + a^2*sigma2
    M_N   = sigma0 + a*M1_N   + a^2*M2_N
    I1_N = (M1_N-sigma1)/sigma0
    I2_N = (M2_N-sigma2)/sigma0
    R_N = 1 + a*I1_N + a^2*(I2_N - (sigma1/sigma0)*I1_N).

On a state without real photons, the measurement z^(N-1) equals one.
The double-virtual term therefore cancels from I2_N. Its shape information
comes from real-virtual and double-real radiation weighted by (F-1), retaining
the proper infrared subtraction between them. The inclusive one-loop rate is
still required in the last normalization term. Dropping that term or integrating
the singular real contributions separately would give an incorrect benchmark.

[McMule's current process table](https://mule-tools.gitlab.io/manual/pieces/)
lists the QED muon-line pieces `ee2mmRFMMMM` and `ee2mmRRMMMM`, together with
`ee2mmFFMMMM`. This makes its subtraction framework a concrete candidate for
the independent check. The initial target should be a pure vector current,
fixed-alpha photonic FSR and identical support, excluding the vacuum-polarization
pieces. This is not a claim that its full Z NNLO interface already exists here:
the table lists the electroweak Z process separately at NLO.

Once that coefficient is available at several masses, fit its L^2, L and
constant pieces and compare directly with the complete leptonff coefficient.
This distinguishes a genuine logarithmic mismatch from a finite remainder
and from the different higher-order resummations. The current independent
alpha scans lack that separation. No new phenomenological theory band or
empirical correction should be inferred from them.

## Reproduction and artifacts

The executable snapshots and campaigns are under `results/z-radiator-audit/`.
Each `kkmc-diagnostic-*` directory contains commands, disjoint seeds, defaults,
client and upstream manifest. The analysis is:

```sh
python3 adapters/kkmcee/diagnose.py \
  --central results/z-radiator-audit/kkmc-central \
            results/z-radiator-audit/kkmc-central-high-stat \
  --coupling results/z-radiator-audit/kkmc-small-alpha \
  --coupling results/z-radiator-audit/kkmc-smaller-alpha \
  --coupling results/z-radiator-audit/kkmc-diagnostic-alpha-half \
  --coupling results/z-radiator-audit/kkmc-diagnostic-alpha-double \
  --small-mass results/z-radiator-audit/kkmc-diagnostic-small-mass \
  --mass-coupling results/z-radiator-audit/kkmc-diagnostic-small-mass-half \
  --mass-coupling results/z-radiator-audit/kkmc-diagnostic-small-mass-double \
  --output results/my-residual-diagnostics
```

`residual-diagnostics/report.json` retains the raw hashes, build identities,
fit covariance, support estimates and results for N=2,3,5,10.
`residual-scaling` and `mass-coupling-scaling` are available as PNG/PDF/SVG.
The plotted coupling panel shows f>=1/2; the weaker alpha/8 and alpha/16
controls remain in the report and all coupling fits. Moments are correlated;
their separate significances or fit p-values must not be combined as independent
tests. Unit tests cover covariance-preserving pooling, fit coefficients and
absolute errors, singular-fit rejection and independent/control-setting guards.
All ten analysis unit tests pass when run directly with the available system
Python. The cached local CTest invocation could not launch: its configured
`build-writing-env/bin/python` points to a removed Python 3.14 installation.
That is an interpreter-path failure, not a successful integrated CTest run.
