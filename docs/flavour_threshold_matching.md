# Charged flavours and time-like threshold matching (v0.21)

The coupled W±/Z match accepts a full charged-fermion basis. The supplied
physical lepton preset contains e, mu and tau, with tau activated near its
mass. The one-family reference and fixed-alpha photonic branch remain
available. This extends the final-state QED calculation; NLO QCD production
still comes from MCFM. Physical hadronic vacuum polarisation and QCD evolution
of coloured fragmentation parents are separate, unfinished inputs.

## Basis, kernels and observations

`QEDFlavourScheme` orders each fermion and antifermion consecutively, followed
by the photon. Each entry specifies its name, charge Q, colour multiplicity
Nc, pole mass and matching scale (zero selects its mass). The observed muon
has Q=-1, Nc=1 and the same mass as `MomentEvolutionParameters::muon_mass`.
The starting scale must include the muon. FF evolution is forward in scale;
the separate coupling transport supports both directions.

Within each active interval, the complete parent matrix evolves the three
observed columns mu−, mu+ and photon. Other charged-parent rows participate
in the evolution and feed back into these columns. The public legacy 3×3
fields select those rows and columns from the full result: their rows do
**not** obey a closed three-species momentum sum in a multi-flavour theory.
Use `flavour_fragmentation` or `flavour_weights()`/`flavour_moment(N)` to
inspect all parents. A full set of observed-flavour columns is not exposed.

Define S2 = sum_active Nc Q² and S4 = sum_active Nc Q⁴. For a=alpha_MSbar/(2π),
da/dln(mu²)=(2 S2/3)a²+S4 a³. The unit-charge time-like distributions are
those of the [existing kernel audit](nlo_photon_kernels.md). Their factors are:

| Channel | Flavour factors at P1 |
| --- | --- |
| fi to fj, pure singlet (each charge) | Nc_j Q_i² Q_j² |
| fi to fi, photonic and identical-charge interference | Q_i⁴ |
| fi to anti-fi, identical-flavour interference | Q_i⁴ |
| fi to fi, nonsinglet pair loop | Q_i² S2 |
| fi to photon | Q_i⁴ |
| photon to fj | Nc_j (Q_j⁴ direct + Q_j² S2 vacuum) |
| photon to photon | S4 |

Produced-flavour colour factors and closed-loop sums are distinct. The
charge/colour decomposition follows the Abelian terms of
[de Florian, Sborlini and Rodrigo, eqs. 51–68](https://arxiv.org/html/1606.02887v3).
Those are space-like kernels; the production distributions use the time-like
continuation documented in the existing audit, rather than substituting the
space-like expressions. Tests check arbitrary charges/colours, per-parent
momentum sums and finite N=1 valence cancellation before evolution.

The boundary generalises the existing O(a) MSbar boundary using each active
flavour's pole mass. Below-threshold parent rows vanish. The target observable
retains the existing inclusive tagged W muon and distinct-parent opposite-charge
Z pair definition. It does not select the leading pair in events with extra
muons, include W decays into other primary lepton flavours, or include tau decays.

## Finite threshold map

For a newly active h at mu_t, let Lh=ln(mu_t²/m_h²), c_h=2 Nc_h Q_h² Lh/3 and
g(z)=[1+(1-z)²]/z. The retained upward map is

```
K_hgamma(z) = Q_h² g(z) [Lh - 1 - 2 ln z]
D_h^above = D_antih^above = a_below K_hgamma ⊗ D_gamma^below
D_gamma^above = (1 - a_below c_h) D_gamma^below
D_light^above = D_light^below
a_above = a_below / (1 - a_below c_h)
```

These are the QED substitutions of the **time-like** NLO relations of
[Cacciari, Nason and Oleari, eqs. 23–27](https://arxiv.org/html/hep-ph/0504192v1).
In particular, the heavy row has a finite term at mu_t=m_h; setting it to
zero would be incorrect. Its analytic Mellin check is

```
K(N)/Q_h² = (N²+N+2)/(N(N²-1)) (Lh-1)
            + 4/(N-1)² - 4/N² + 2/(N+1)² .
```

The rational completion of coupling matching has the same one-loop expansion
and makes a D_gamma continuous for this map. It does not supply unknown
two-loop matching constants. Downward coupling matching uses its exact inverse.
Simultaneous thresholds use the same below-threshold photon row and the sum
of their photon/coupling logarithms. Explicit matching scales are restricted
to [m/4,4m]; this is a perturbative matching prescription near each mass.

## Input scheme and hard matching

`alpha_msbar_at_mu0` in the C++ API is always an MSbar input. The optional
`leptonic_alpha_msbar` conversion uses

```
1/alpha_MSbar(mu) = 1/alpha_OS
                   - sum_active Nc Q² ln(mu²/m²)/(3π) .
```

It is the rational one-loop conversion; higher conversion coefficients are
not included. The beta-function and OS/MSbar conventions are consistent with
[Baikov et al.](https://arxiv.org/html/1207.2199v2). Automatic conversion rejects
any coloured entry: perturbative light-quark pole-mass logarithms are not used
as a substitute for hadronic vacuum polarisation.

The adapter's `--flavours leptons` defaults to on-shell input; legacy and
custom profiles default to explicit MSbar input. The numeric default alpha
is the existing reference 1/137.035999084. The preset masses in GeV are
e=0.00051099895069 ([CODATA 2022](https://www.physics.nist.gov/cuu/pdf/wall_2022.pdf)),
mu=0.1056583755 (retained reference), tau=1.77693
([PDG 2025](https://pdg.lbl.gov/2025/download/db2024.pdf)). Threshold variation
moves only tau in the built-in preset. Each run logs the actual masses,
matching scales, converted alpha(mu0), alpha(mu_frag) and running logarithm.

The finite hard decay uses the same threshold-matched alpha at the original
hard invariant mass M. Its retained C2 expansion replaces the one-family
beta0 logarithm by B(M)=integral beta0 dln(mu²)+sum_crossed c_h. C1/C2 include
the same finite threshold maps and the same lattice measure as the resummation.
Coefficients expand in a_MSbar(mu0), even when the adapter accepts OS input.
An expansion in a_OS would also transform C2 by the one-loop input-conversion
coefficient times C1. No two-loop finite boundary/threshold accuracy is claimed.
Fixed coupling is a mathematical control; two-loop running is the physical
default of the new mode.

## Running and verification

```sh
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-flavours-v021/mcfm_leptonff --output results/my_wplus \
  --coupled-hard --fraction-cells --hard-samples 4 --spectra --threads 1 \
  --flavours leptons
```

Use `--threshold-factor .5` or `2`, or `--alpha-scheme msbar` for an explicit
MSbar input at mu0. A custom file has strict lines
`name charge colours pole_mass_GeV matching_scale_GeV`, comments starting `#`,
and an observed entry named `mu`. See [lepton input](../configs/leptons_thresholds.dat).
`--flavours custom --flavour-file FILE --alpha-scheme on-shell` accepts only
uncoloured leptonic inputs. The
[charged-fermion test file](../configs/charged_fermion_threshold_test.dat)
is a synthetic pure-QED charge/colour control, not a Standard Model prediction.

The [v0.21 benchmark](../benchmarks/flavour_threshold_reference.md) records
the independent solver, core checks, full host controls and residual warnings.
PHOTOS comparisons with this mode remain diagnostic because PHOTOS uses a
fixed coupling without pair channels. The aligned photonic branch is unchanged.
