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

`adapters/dyturbo/` contains original MIT-licensed client/build code. Its optional
builder downloads checksum-pinned upstream archives into an ignored directory,
retains upstream licence files, and makes documented build compatibility edits
to isolated HELLx/DYTurbo sources. DYTurbo is GPLv3; linked native executables
remain subject to the backend and dependency licence terms. These executables,
downloaded sources and PDF grids are not distributed with the MIT library.

`adapters/renesance/` contains original MIT-licensed reference-driver code.
Its builder obtains the checksum-pinned GPLv3 ReneSANCe 1.3.11 release and
keeps its SANC modules, generated FSR-only integrated expressions, and linked
executables in an ignored build directory. The generated expressions retain
the SANC coefficients and copyright; they and the linked reference executable
are subject to the external licence terms. They are not part of the MIT source
distribution. PYTHIA and PHOTOS remain separately obtained dependencies.

`adapters/kkmcee/` contains original MIT-licensed build/reference code for
the separately obtained official KKMCee v4.32.01 Fortran release. The builder
pins its Git commit, retains the upstream files and builds the original
amplitudes without copying them into leptonff. KKMCee and its bundled
TAUOLA/PHOTOS/utility components retain their respective copyright, licensing
and citation terms. Their sources and linked executables stay in ignored
external/build/result directories and are not relicensed or distributed here.

`adapters/mcmule/` contains original MIT-licensed API-client, analysis and
importance-coordinate code for separately obtained McMule v0.8.0. Its builder
pins the commit and recursive submodules. A documented preparation script adds
an explicit call to the original two-body phase-space routine to avoid an
indirect-call compiler/platform failure; it edits no amplitude or subtraction
formula. Prepared upstream sources retain GPL-3.0-or-later terms in the ignored
external directory. McMule, OpenLoops, COLLIER and the other dependencies retain
their own licence and citation requirements. Their sources and linked reference
executables are not distributed as part of the standalone MIT library.

MCFM, PHOTOS, PYTHIA, KKMCee, McMule, ReneSANCe/SANC, RECOLA, LHAPDF/PDF sets and DYTurbo are optional external
software or data. Their sources, libraries and datasets are not bundled.
They retain their own licensing and citation requirements; the root MIT
licence does not relicense them or executables combining them.
The standalone library and core tests require none of these packages.

The example input `configs/hadronic_alpha_pdg2025.dat` records its PDG source
and convention. Physics references and implementation conventions are in
`docs/`; the software citation does not replace those scientific citations.
