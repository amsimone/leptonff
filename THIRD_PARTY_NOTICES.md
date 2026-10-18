# Licensing and external dependencies

The root MIT licence covers original leptonff code and documentation except
where a file or directory states otherwise. Copyright: 2026 Simone Amoroso.

## MCFM adapter

`adapters/mcfm/` is distributed under GPL-3.0-or-later; see its
[notice](adapters/mcfm/LICENSE) and the
[complete GPL text](LICENSES/GPL-3.0-or-later.txt).
The replacement `nplotter_user.f90` adapts the MCFM 10.3.1 user hook:
Copyright (C) 2019–2022, respective authors of MCFM. The relinker also
modifies local MCFM integration sources when building an isolated host
executable. Generated host sources and linked executables are not part of
the standalone MIT library or this source distribution.

## Separately obtained software

MCFM, PHOTOS, RECOLA, LHAPDF/PDF sets and DYTurbo are optional external
software or data. Their sources, libraries and datasets are not bundled.
They retain their own licensing and citation requirements; the root MIT
licence does not relicense them or executables combining them.
The standalone library and core tests require none of these packages.

The example input `configs/hadronic_alpha_pdg2025.dat` records its PDG source
and convention. Physics references and implementation conventions are in
`docs/`; the software citation does not replace those scientific citations.
