# Research extensions: fixed order, endpoint and dressed leptons

Research notes, 18 October 2026. These are proposed extensions, not implemented
accuracy claims. The present matched photonic prediction is NLO in the decay
and NLL in lepton-mass logarithms. The McMule study independently tests our
second-order logarithmic coefficients for normalized inclusive vector-current
moments; it has not upgraded the production prediction to NNLO.

The objective throughout is to validate and improve **leptonff**. Independent
codes provide reference coefficients and measurements over their common scope;
explaining every difference between external generators is not a prerequisite.
The following proposals separate established ingredients, derivations still
needed, and implementation gates. Factorization formulas labelled schematic
are organizing proposals rather than results already derived for this package.

## Four distinct extensions

| Extension | New information | First useful target | Present status |
|---|---|---|---|
| NNLO QED + collinear NLL | Complete finite second-order contribution for a specified decay observable | Angle-inclusive neutral-current bare pair mass | Independent moments exist; differential matching is missing |
| Endpoint resummation | Logarithms enhanced as the specified energy loss tends to zero | Fixed-Q bare pair-mass cumulatives | Part of the soft–collinear tower is already in DGLAP; uniform endpoint matching is missing |
| Soft plus collinear resummation | Coherent wide-angle radiation with its collinear overlap removed | On-shell Z decay with a global soft-energy measurement | Existing hard currents and FFs are reusable; a measurement-dependent soft operator is missing |
| Dressed-electron cone plus veto resummation | Cone, veto and mixed logarithms, with a physical veto efficiency | Two dressed electrons and total photon energy outside the cones | Needs a new angular measurement and factorization contract |

These are related but not interchangeable. Completing a constant at second
order does not resum endpoint logarithms. Adding soft exponentiation does not
automatically establish collinear NLL. Dressing does not mean replacing the
muon mass by the electron mass in a bare-lepton fragmentation function.

## 1. NNLO matching: what transfers between Z and W

In a schematic single-flavour notation the decay contains

$$
d\Gamma_Z=H_Z\otimes D_-\otimes D_+,
\qquad d\Gamma_W=H_W\otimes D_\ell.
$$

The fragmentation functions are universal. The hard operators, including
their angular/spin dependence and finite terms, depend on the process.
At second order the full answer contains the two-loop fragmentation boundary,
the second-order hard coefficient, and products of lower-order terms. Those
products are not all retained by our current single-subleading-insertion NLL
prescription. They must be present in a complete fixed-order overlap.

For a common normalized observable, additive matching would read

$$
R_{\mathrm{NNLO+NLL_c}}=R_{\mathrm{NLL_c}}
+R_{\mathrm{FO},\leq\alpha^2}
-[R_{\mathrm{NLL_c}}]_{\leq\alpha^2}.
$$

Here the subscript c means collinear mass-log accuracy. Couplings, support,
branching normalization, mass treatment and measurement must agree. The same
scheme must be expanded in the subtraction. A massive fixed-order input can
restore its mass powers through the matched order without claiming to resum
those powers. Additive matching of fixed order and resummation has a close
fragmentation precedent in [Bonino, Cacciari and Stagnitto, section 2](https://arxiv.org/html/2312.12519).

### Neutral current

Our McMule extraction gives the second-order **normalized vector-current
shape**, not a universal two-loop FF. Its finite remainder mixes hard,
boundary and lower-order-product contributions. Moments alone cannot separate
those contributions or reconstruct the distribution uniquely.

The concrete first extension is a conservative differential/cumulative
response in the bare pair-mass fraction. It needs a controlled endpoint
representation, covariance and a precise overlap subtraction. For the
angle-inclusive, parity-even Z/gamma* application, establish explicitly the
vector/axial correspondence in the no-fermion-loop leading-power sector.
Do not transfer exact vector mass-power terms to an axial current without
checking that mapping. Fiducial cuts and angular observables additionally need
the corresponding differential decay response.

For the normalized fixed-Q mass measurement F(Born)=1, double-virtual terms
cancel in the F-1 action. The inclusive NLO normalization counterterm still
enters the second-order ratio. An absolute NNLO partial width requires the
missing two-loop normalization. The cancellation cannot be assumed for a new
fiducial measurement: derive it using that measurement and spin structure.

### Charged current

A W NNLO match needs a dedicated second-order W decay calculation. The charged
W participates in radiation and interference; its hard correction is not half
the neutral-current correction. Universal FF improvements transfer, while the
finite W hard coefficient does not come from the present McMule setup.
Our existing gauge-consistent first-order W current is a starting point for
matching and soft-limit checks, not a substitute for that calculation.

This is photonic QED accuracy **in the decay**. It is separate from DYTurbo's
NNLO+NNLL QCD production accuracy and from full NNLO electroweak or mixed
QCD–electroweak corrections. A full NNLL collinear extension would also need
the requisite higher evolution, coupling and matching ingredients; NNLO
fixed-order matching alone does not provide it.

## 2. Define which endpoint is being resummed

Use distinct symbols rather than calling every fraction z:

| Variable | Meaning | Endpoint |
|---|---|---|
| $x_i$ | Fragmentation fraction on charged leg i | $x_i\to1$ |
| $z_m=m_{\ell\ell}^2/Q^2$ | Bare pair mass at fixed parent Q | All radiated energy becomes soft |
| $z_d=m_{\mathrm{dressed}\,\ell\ell}^2/Q^2$ | Pair mass after photon recombination | Energy outside the dressing cones becomes soft |
| $\beta=E_v/Q$ | Explicit energy-veto fraction | $\beta\ll1$ |

A lepton pT Jacobian edge is not identical to any of these limits. Its
prediction also depends on production recoil, decay angle and fiducial cuts;
QED endpoint resummation does not replace QCD recoil resummation.

Our current DGLAP exponent already resums soft–collinear endpoint-enhanced
terms: the leading kernel has a whole-plus $1/(1-x)$ singularity. The endpoint
atom and soft cutoff in the numerical FFT representation are regulators of
that calculation, not new physical resummation scales or a physical photon
veto. The older endpoint-panel work improves quadrature, not logarithmic
accuracy.

### A diagnostic derived from our present formula

For the inclusive mass observable at mu_frag=Q, define a=alpha/(2*pi),
L=log(Q^2/m^2), and Mellin N conjugate to z_m. Our first-order coefficient is

$$
r_{1,N}=2(LP_N+C_N),\quad
P_N=\int_0^1 dx\,p(x)(x^{N-1}-1),\quad
C_N=\int_0^1 dx\,p(x)(\ln x-1)(x^{N-1}-1),
$$

with p(x)=(1+x^2)/(1-x). For positive integer N,

$$
P_N=\frac32-2H_N+\frac1{N(N+1)},\qquad
C_N+P_N=2H_{N-1}^{(2)}+\frac1{N^2}+\frac1{(N+1)^2}-\frac54.
$$

Consequently, for large N and $\bar N=N e^{\gamma_E}$,

$$
P_N=-2\ln\bar N+\frac32+O(1/N),\qquad
C_N=-P_N+\frac{\pi^2}{3}-\frac54+O(1/N),
$$

and the coefficient of ln(N-bar) in r1 is -4*(L-1). The double endpoint
logarithms of the separate MSbar boundary and hard coefficient have canceled
in this **complete observable** at first order.

The current collinear result contains an LL factor with large-N behaviour
exp[-4*a*L*ln(N-bar)] and a single finite insertion supplying +4*a*ln(N-bar)
at first order. Exponentiating the physical soft radiation would instead
organize that leading soft exponent with L-1. The higher powers generated by
the finite soft contribution are beyond the present single-insertion
collinear counting but can matter when ln(N) itself is large.

This derivation identifies a useful extension. It is not a derivation of the
complete endpoint exponent: angular soft interference, constants, mass powers
and the measurement must also be treated consistently. The existing
`quadratic_finite` variation includes one such higher-order term, but is not
an endpoint-resummed prediction.

### Proposed endpoint calculation

Derive the physical large-N radiator, including its hard coefficient and
boundary, and match it to the current collinear result with their common
soft–collinear terms removed. Preserve the complete first-order coefficient;
where available, preserve the independent second-order coefficient as well.
Use a profile to turn off the endpoint approximation away from the endpoint.
Check normalization of the matched full-support prediction rather than
extending a large-N formula to N=1 without justification.

[Cacciari–Catani](https://arxiv.org/abs/hep-ph/0107138) supplies the classic
fragmentation/threshold organization. Its QCD ingredients must be reduced to
the specified QED sector and rederived for our observable, rather than copied
by replacing a coupling. [Ghira–Marzani–Ridolfi](https://arxiv.org/html/2309.06139)
shows why soft and small-mass limits require care when combined. Audit the
regions where endpoint scales cross m or m^2/Q; do not assume a fixed-z
leading-power expansion is uniform in every endpoint hierarchy.

The electron-PDF large-z analyses of
[Bertone et al.](https://arxiv.org/html/1911.12040) and
[Frixione](https://arxiv.org/html/2105.06688) are useful for analytic/numerical
matching and scheme sensitivity. They concern spacelike PDFs. A timelike FF
extension requires transformed kernels, boundaries and hard coefficients
together; the endpoint shape of an unphysical FF is not an observable test.

## 3. Coherent soft radiation plus collinear evolution

Collinear fragmentation is universal per hard direction. Wide-angle soft
radiation depends on the full set of charged directions through interference.
For an on-shell neutral decay there are two outgoing charged legs; for a W
decay there is an incoming charged W and an outgoing charged lepton.

The leading soft current has the form

$$
J^\mu(k)=\sum_i\eta_i q_i\frac{p_i^\mu}{p_i\cdot k},\qquad
\sum_i\eta_i q_i=0,
$$

where eta is -1 for incoming and +1 for outgoing charges. Its squared current,
including interference and massive self terms where appropriate, defines
the emission measure. Independent leg probabilities alone do not supply this
angular information. Soft exponentiation is established in QED; the
[CEEX formulation](https://arxiv.org/abs/hep-ph/0006359) is a relevant
amplitude-based reference, not a replacement for our collinear bookkeeping.

A proposed scalar benchmark organization is

$$
\ln\widetilde R_{\mathrm{joint}}
=\ln\widetilde R_{\mathrm{collinear}}
+\ln\widetilde R_{\mathrm{soft}}
-\ln\widetilde R_{\mathrm{overlap}}
+\text{matching terms}.
$$

The transforms and overlap must be derived for the same measurement. This
equation is schematic and is not a prescription for taking scalar logarithms
of the full flavour-matrix evolution. Multiplying the existing FF by a YFS
factor without the overlap subtraction would count soft–collinear radiation
twice. In effective-theory language the overlap is removed by the appropriate
zero-bin/subtracted soft and collinear functions.

An elementary check is the real-kernel decomposition

$$p(x)=\frac{2}{1-x}-(1+x).$$

The first term is already present in the DGLAP evolution. It cannot be added
again as a new soft probability. The regular remainder in this algebraic
decomposition is signed; it is not by itself a positive branching rate.

For an additive measurement and fixed charged directions, independent eikonal
photons in fixed-alpha QED without fermion loops have a particularly simple
generating functional:

$$
\ln\widetilde S(\nu)=\int dP_{\mathrm{eik}}(k)
\left[e^{-\nu v(k)}-1\right].
$$

The virtual term gives S-tilde(0)=1. Multiple-emission inversion is still
required: a veto on the **sum** of photon energies is not a veto on each
individual photon energy. Recoil, hard emissions, running/pair corrections
and more complicated clustering go beyond this restricted independent-soft
statement. [Dai–Kim–Leibovich](https://arxiv.org/html/2103.03963) gives a QED
example separating hard, soft, massive collinear and collinear-soft scales in
exclusive leptonic decays; it is a useful starting point for the bare veto
benchmark, not an existing derivation for every W/Z fiducial measurement.

Keep accuracy labels separate: fixed-order QED, collinear logarithms in L,
endpoint/veto logarithms, and cone logarithms. State which terms are retained
in the exponent and which mixed powers are counted as leading. An unqualified
second “NLL” would otherwise hide a change of expansion.

## 4. A precise dressed-electron plus veto benchmark

Begin in a fixed-Q on-shell neutral-current decay, fixed alpha, photons only,
with no ISR, production shower, pairs or detector response. Define two
nonoverlapping cones in the parent rest frame and use four-vector recombination,

$$
p_{e,i}^{\mathrm{dress}}=p_{e,i}
+\sum_{\gamma\in C_i(R)}k_\gamma,\qquad
E_{\mathrm{out}}=\sum_{\gamma\notin C_1\cup C_2} k_\gamma^0.
$$

Accept events with E_out<E_v, inclusive in all photon energies inside the
cones. The first analytic benchmark may fix the axes to projected Born
directions, with that projection explicitly part of the definition. A final-
bare-electron axis or a recoil-free clustered axis is a different definition
and needs its own matching. Do not update cone axes implicitly while summing
photons. Specify assignment of photons if cones can overlap.

The physical hadron-collider extension would use lab delta-R cones, a declared
axis, transverse rather than rest-frame energy if desired, and finite detector
acceptance. A rest-frame angle is not automatically a lab delta-R. The angular
map and veto frame must be retained in the production response.

An annular isolation requirement, a hardest-photon veto, and a sum-energy veto
are distinct observables. The latter two agree for one photon but differ
beginning with two emissions. A hadronic isolation veto additionally measures
QCD radiation and cannot be supplied by a decay-only QED soft function.

### Which logarithms survive dressing?

For fixed R well above m_e/E and an observable inclusive in collinear
electron–photon recombination, the leading electron-mass singularities cancel.
The radiation lost **outside** the cone contains logarithms of R instead.
Thus a pure photon-sector dressed observable is not predicted by merely
evolving a bare electron FF down to the electron mass. A track-energy
requirement or a separate bare-electron cut can spoil that inclusiveness and
restore mass-log sensitivity.

For a narrow cone and a tight veto the one-emission loss contains the region

$$
\int_{E_v}^{E}\frac{d\omega}{\omega}
\int_{R^2}^{O(1)}\frac{d\theta^2}{\theta^2}
\ \sim\ \ln\frac{E}{E_v}\ln\frac1{R^2}.
$$

This is a logarithmic-region estimate, not the complete coefficient: charges,
two-leg interference and finite angular geometry must be included. An
exclusive efficiency is Sudakov suppressed as the veto tightens. Independent
evolution of two unrelated logarithms is not a sufficient joint calculation.
For a concrete QED treatment of cone-inclusive real radiation and the
small-mass cancellation, see
[Tomalak et al., sections 4.4 and 6.2](https://arxiv.org/html/2204.11379).
Its neutrino hard process differs from ours, so transfer the universal
functions only after checking conventions and the observable.

### Candidate factorization and scales

For the initial small-cone/soft-veto hierarchy, a useful organizing proposal is

$$
d\Gamma[F]\sim H_{ab}
\left[\prod_i J_{e,i}\otimes S_{\mathrm{cs},i}\right]
\otimes S_{\mathrm{wide}}[F].
$$

The measured electron-jet functions describe cone-inclusive collinear
radiation; the collinear-soft functions resolve soft radiation near a cone
boundary; the wide-angle function includes the charged-direction interference
and the veto geometry. All functions need compatible overlap subtraction.
They need not depend on independent variables, and some may be combined in
another equally valid factorization. This is a proposed decomposition to
derive, not a theorem already demonstrated for our dressing algorithm.

For E~Q, m_e much smaller than E*R, and E_v/Q much smaller than R much smaller
than one, characteristic virtuality scales include

$$
\mu_H\sim Q,\quad \mu_J\sim E R,\quad
\mu_S\sim E_v,\quad\mu_{\mathrm{cs}}\sim E_v R.
$$

The choice E_v much smaller than E*R is an initial regime, not a requirement
for every physical event. Other hierarchies need compatible profiles or
matching. In the bare massive veto problem the analogous collinear and
collinear-soft scales instead involve m and E_v*m/E. The R~m_e/E transition
therefore needs a massive treatment; the massless cone formula cannot be
continued to R=0 to recover a bare-electron result.

Check cancellation of the summed anomalous dimensions, with all collinear-
soft pieces included. Running alpha, flavour thresholds and vacuum
polarization must be added consistently when departing from the fixed-alpha
photonic benchmark. A QCD infrared Landau-pole prescription is not a QED
ingredient to import automatically.

### Endpoint and veto can measure the same radiation

Let K_out be the total momentum outside both cones. Four-vector dressing gives

$$
z_d=\frac{(P-K_{\mathrm{out}})^2}{Q^2}
=1-\frac{2E_{\mathrm{out}}}{Q}+\frac{K_{\mathrm{out}}^2}{Q^2}.
$$

Consequently, close to the endpoint, 1-z_d and E_out are the same measurement
at leading soft power. Their cuts are correlated. Multiplying an endpoint
Sudakov by a second independent outside-energy Sudakov would double count.
A joint transform must encode both measurements on each photon; here their
leading-power dependence reduces to one combination of transform variables.
Different single-electron observables or an annular veto can supply genuinely
additional measurement dependence.

### A veto efficiency must not be normalized away

Our current inclusive decay response preserves a supplied branching
normalization. With a veto, the integrated response should instead return

$$\epsilon(R,E_v)=\Gamma[E_{\mathrm{out}}<E_v]/\Gamma_{\mathrm{incl}}.$$

Keep this efficiency and the conditional accepted-event shape separately.
Forcing the vetoed response to integrate to one would erase a physical effect.
The no-veto limit is recovered by a full matched prediction with resummation
profiles; extrapolating a small-E_v asymptotic formula to E_v~Q is not that
check.

## 5. Angular restrictions, clustering and pair effects

An outside-cone veto is geometrically non-global, but one should not attach a
standard QCD non-global evolution factor without examining the theory. In the
fixed-alpha, no-fermion-loop, independent-photon eikonal limit with fixed axes,
photons do not radiate softer photons: the familiar non-Abelian secondary-
emission mechanism is absent. Geometry still enters the one-photon integral,
and a sum-energy cut still needs the correct multi-emission inversion.

This simplification is restricted. Recoil or a clustering algorithm can
correlate the measurement even in an Abelian theory. The QCD analysis of
[Kelley–Walsh–Zuberi](https://arxiv.org/abs/1202.2361) demonstrates Abelian
clustering logarithms; it motivates a two-photon check for the chosen QED
algorithm, not a claim of an identical coefficient. Fermion-pair production
introduces additional charged emitters and a separate physical-content issue.

For electrons, a prediction including gamma -> e+e- must specify how conversion
products, charge and electron identity enter dressing and vetoes. Photon-only
dressing need not have the same collinear inclusiveness once a photon is
resolved into a pair. Either keep the relevant masses with an explicit
particle-level selection, or derive a suitable flavour-safe clustering
definition. Simply turning on an inclusive pair FF does not define an
exclusive electron event. [Flower–Schönherr](https://arxiv.org/html/2210.07007)
studies this sensitivity in Z -> e+e- and compares dressing strategies.
Its quoted effects are not numerical predictions for our proposed veto.

At hadron level, initial-state charged partons also contribute to the soft
current. A decay-only calculation deliberately excludes their radiation and
initial–final interference. This restriction is especially important for a
charged resonance: if the veto probes energies comparable to its width, the
production–decay separation and nonfactorizable radiation require assessment.
An on-shell decay benchmark is well-defined, but does not establish a complete
vetoed pp prediction in that region.

## 6. Implementation implications and validation gates

Reuse the existing distributions, finite hard currents, spin tensors,
production records, covariance and matching-expansion machinery. Add new
objects alongside the inclusive FF rather than changing its interpretation:

| Proposed object | Information it must retain |
|---|---|
| Dressing/veto measurement | Cone metric, axes and frame; recombination; sum/max veto; acceptance; pair prescription |
| Electron-jet response | Cone and mass dependence, energy fraction, and recoil information required by the measurement |
| Soft response | Charged directions, charges and incoming/outgoing signs; measured energy; interference |
| Matched joint response | Explicit overlap and fixed-order expansion; separate inclusive normalization and veto efficiency |
| Response cache/archive metadata | Every physical measurement parameter and the soft geometry; numerical profiles and precision separately |

A bare longitudinal FF response table cannot reconstruct arbitrary photon
angles after the fact. Existing production archives may be reusable if they
retain the needed parent kinematics and spin information, but new QED response
tables are required. A hadronic veto additionally needs production radiation
information, beyond a record containing only parent production tensors.

The initial validation sequence should be:

1. **Analytic endpoint audit.** Verify the displayed large-N identities,
   cancellation of boundary/hard endpoint logs and the physical first-order
   soft-current coefficient. Distinguish mass-log and endpoint counting.
2. **Exact first-order dressed response.** Integrate the existing hard current
   with the chosen photon-angle/recombination/veto measurement, including
   subtraction counterevents. Check soft and collinear cancellation, Ward
   identities and independence of slicing/subtraction parameters.
3. **Factorization/overlap check.** Expand the proposed soft, jet and
   collinear-soft product through first order and reproduce that response.
   Check all mixed logs and scale cancellation. A plausible Sudakov shape is
   not a substitute for this equality.
4. **Independent second-order tests.** At fixed Q, compare moments/cumulatives
   and veto efficiencies with the corresponding McMule measurement where
   supported, targeting our claimed logarithmic coefficients. Check two-
   photon sum versus max vetoes and clustering explicitly. More endpoint
   statistics will be needed than in the current moment study.
5. **Limit and numerical scans.** Vary R, E_v, m_e, soft and grid cutoffs,
   transform inversion and profiles independently. In the photon-only
   collinear-safe regime, mass sensitivity should disappear at leading power.
   Recover bare limits with a massive matched treatment and no-veto limits
   with full fixed-order matching, not outside a formula's valid hierarchy.
6. **Boosted and charged-current extension.** Test lab cone geometry and
   production recoil first, then W charge/spin and charged-parent soft
   emission. Use a dedicated W fixed-order anchor. Full NNLO W matching is a
   separate gate from adding soft/veto logarithms to the existing NLO W match.

For plots, separate inclusive rates, veto efficiencies and normalized accepted
shapes. Useful targets are fixed-Q bare/dressed pair-mass cumulatives,
efficiency versus E_v at several radii, and an (R,E_v) map of corrections.
Only after these pass should fiducial electron pT/eta spectra be folded with
QCD production. Quote statistical, discretization, inversion/profile and
perturbative-scale effects separately.

## 7. Suggested order of work

Soft/collinear factorization, QED exponentiation and threshold-improved
fragmentation already have substantial literature. The potential contribution
here is a consistent combination for the declared W/Z measurements, explicit
fixed-order validation, and a reusable spin-dependent production interface.
A new factorization or novelty claim would require a separate derivation and
a closer comparison with that literature.

First complete the inclusive neutral-current NNLO matching response using the
existing McMule infrastructure, while deriving and auditing the endpoint
structure of our current formula. These are useful independent deliverables.
Next build a first-order dressed-electron cone-plus-total-energy-veto
reference and derive its overlap-subtracted soft/collinear factorization.
Implement a joint resummation only after that expansion closes. Start with
fixed-alpha photons and a clean on-shell neutral decay; then add realistic
lab cones, W decay, running/pairs and production–decay effects in separately
labelled steps.

The dressed-veto branch is scientifically attractive because it extends the
set of predictable observables, rather than only adding a small finite term
to an existing moment. Its numerical importance should first be measured by
an NLO (R,E_v) scan; the electron's smaller mass alone does not predict that
importance after inclusive dressing. A simple fixed-cone Abelian soft benchmark
may be evaluated with angular quadrature and transform inversion, then cached
per charged geometry. That could avoid a large exclusive multi-photon MC
campaign for each study, but its speed and precision have not been benchmarked.

This program can improve Z earlier than W at complete NNLO fixed order, while
shared soft/jet components can benefit both at their established matching
order. None of these notes changes the accuracy label of the current code.
