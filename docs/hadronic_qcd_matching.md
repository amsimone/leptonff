# Hadronic input and coloured-parent evolution (v0.22)

The `hadronic-pointlike` profile adds a sourced hadronic alpha input and LO
QCD quark/gluon evolution to the NLO QCD W±/Z fragmentation benchmark. Its
light-quark fragmentation boundary is an **explicit pointlike approximation**.
This is an optional combined-evolution benchmark with a sourced
vacuum-polarisation input. Mixed kernels and a fitted boundary concern that
larger programme; they are not prerequisites to the intended pure-QED NLL
calculation on NLO QCD production. `--qcd-off` disables strong FF evolution
while retaining the explicit pointlike/HVP prescription. See the
[accuracy audit](qed_accuracy_audit.md) for the scoped comparison and its limits.

## Input and coupling matching

[The PDG 2025 standard-model review, section 10.3.1](https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf)
uses the timelike, on-shell, five-flavour input
Delta alpha_had5(MZ²)=0.02760 ± 0.00010 with alpha_s(MZ)=0.1193.
The checked-in [input file](../configs/hadronic_alpha_pdg2025.dat) records these
values and MZ=91.1876 GeV. This Delta is neither an MSbar coupling nor a real
hadronic-pair fragmentation distribution. Its uncertainty is one input
contribution, not a complete prediction uncertainty.

For alpha0=alpha_OS and a_s=alpha_s/(2π), the high-scale fermionic QED anchor is

```
1/alpha_MSbar,fermionic(MZ) = 1/alpha_MSbar,leptonic(MZ)
  - Delta alpha_had5,OS(MZ²)/alpha0
  - (1/pi) [55/27 + (alpha_s(MZ)/pi)(605/108 - (44/9) zeta(3))] .
```

The leptonic term uses the existing one-loop OS input at mu0=m_mu and
two-loop, threshold-matched leptonic running. The hadronic finite conversion
uses the quark part of PDG eq. (10.21): 55/27 is the five-quark contribution
to its fermionic 100/27. **W and top loops are excluded**. Higher QCD conversion
terms, higher finite QED input matching and quark mass power terms are omitted.
The result at the default anchor is alpha=0.00782681977749425
(inverse 127.76581401241222). It should not be identified with PDG's full
electroweak MSbar coupling, which includes additional contributions.

An extra matching coefficient c_H is fitted at mu_H to reproduce this anchor
after the actual perturbative quark running and ordinary heavy thresholds.
It is a residual coefficient in this truncated scheme, rather than a measured
local polarisation function. At the matching scale,

```
a_above = a_below / (1 - a_below c_H)
D_gamma,above = (1 - a_below c_H) D_gamma,below .
```

Ordinary heavy-threshold photon coefficients at the same scale add to c_H.
The finite photon map, hard-scale coupling, running logarithm and C1/C2
expansion all use the same coefficient. Fitting the residual instead of
adding the full measured Delta on top of perturbative quark running avoids
double counting. When varying alpha for a perturbative expansion, hold c_H
fixed after calibration. The zero-alpha control uses the physical reference
normalisation to define c_H and then sets the actual input coupling to zero.

The installed Jegerlehner `hadr5x19.f` routine supplies an independent alternative
OS HVP anchor. [Its documentation](https://people.physik.hu-berlin.de/~fjeger/alphaQED19.pdf)
describes the spacelike/timelike conventions. `python/check_external_hvp.py`
compiles that external file read-only, hashes it and exports its output as a
separate input. It does not vendor or modify the installed implementation.
The double-argument wrapper still uses single-precision interpolation internally.

## Partons, boundary and QCD accuracy

The parent basis adds a gluon after the photon. Existing fermion indices and
the selected mu−/mu+/photon outputs remain stable. All parent rows evolve;
the selected 3×3 output is not a closed momentum system. The default profile
contains e/mu/tau and u/d/s/c/b with their physical charges and Nc=3 for quarks.

The u/d/s masses are zero and their matching scale is mu_H. Their
fragmentation boundary into the selected observed species is zero at that
scale. No light-quark pole mass or finite heavy-particle matching formula is
assigned to them. The `pointlike` flag is explicit in the input type and logs.
For a custom flavour file, an optional sixth field `pointlike` selects this
massless boundary. It is permitted only for Nc=3 and an explicit matching
scale of at least 1 GeV.

Charm and bottom use the explicit reference pole-mass choices 1.67 and
4.78 GeV; charm is activated at max(mu_H,1.67 GeV). Their QED threshold map is
the finite time-like map of v0.21. The default mu_H is 2 GeV; variation between
1 and 4 GeV is permitted. There is no perturbative QCD evolution below mu_H.
The zero pointlike boundary is a model choice whose replacement by physical
nonperturbative matching information remains open. It does not include tau
decays, hadron-decay muons or additional primary W decay flavours.

The time-like LO QCD matrix in b=alpha_s/(2π) contains
q→q with CF Pff, q→g with CF Pfg, g→q/anti-q with TR Pgf, and the usual
gluon self kernel, with CF=4/3, TR=1/2 and CA=3. Per-parent N=2 momentum sums
vanish. The coupling obeys

```
db/dln(mu²) = -(11/2 - nf/3) b² .
```

It runs piecewise from the reference input and is continuous across QCD
thresholds at LO. The QCD heavy-parent boundary is zero at this accuracy;
no NLO QCD threshold or finite QCD boundary accuracy is claimed. The independent
time-like reference of [Bertone et al.](https://arxiv.org/html/1501.00494v2)
provides the convention and evolution framework.

The implemented generator is b P_QCD0 + a P_QED0 + a² P_QED1.
Mixed a b splitting kernels and mixed beta terms are **absent**. Iterating
this generator resums QCD dressing of the QED parents but does not supply
those missing mixed coefficients. Production alpha_s/PDF choices are recorded
separately from the fragmentation reference input. They need alignment for
precision phenomenology.

C1/C2 are the expansion in a_MSbar(mu0) with the strong coupling held fixed
as an independent input. They satisfy the corresponding QCD-dressed
inhomogeneous ODEs and include the finite matching maps. Strong evolution is
not approximated by its first insertion in these coefficients. The old
analytic polynomial integration remains in use when QCD is absent.

## Use

```sh
python3 adapters/mcfm/run_reference.py /path/to/MCFM-10.3.1 wplus \
  --adapter build-mcfm-hadronic-final-v022-release/mcfm_leptonff \
  --output results/my_hadronic_wplus --threads 1 --spectra \
  --coupled-hard --fraction-cells --hard-samples 4 \
  --flavours hadronic-pointlike
```

This defaults to the sourced input above, on-shell alpha0 and a lattice log
span of 64. Other modes retain their previous defaults and remain regression
references. Available controls are `--hadronic-input FILE`,
`--hadronic-shift -1`/`1` (input standard deviations), `--hadronic-scale SCALE`,
`--alpha-s VALUE`, `--qcd-off` and `--log-span VALUE`. `--qcd-off` disables the
QCD fragmentation evolution while preserving the same hadronic coupling
anchor and conversion input. Input files contain all four fields
`reference_scale`, `delta_had5_os`, `uncertainty`, `alpha_s_reference`.
All are validated, copied into each run and hashed.

The [benchmark](../benchmarks/hadronic_qcd_reference.md) records the independent
checks, host campaign and limitations. PHOTOS comparisons remain diagnostic.
