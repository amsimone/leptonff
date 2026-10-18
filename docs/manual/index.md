# leptonff scientific manual

leptonff computes perturbative timelike QED fragmentation and matched
final-state radiation for bare-lepton W/Z measurements. The reusable objects
are distributions, evolution operators and spin-correlated decay responses.

The C++17 core runs independently of a QCD generator. Matched W/Z modules add
finite hard-decay terms. Generator adapters supply QCD production; production
archives and replay support repeated studies of compatible QED models and cuts.

```{toctree}
:maxdepth: 2

installation
physics
tutorials
interfaces
validation
architecture
python
```

Start with [installation](installation.md) and the [three tutorials](tutorials.md).
Read the [scientific contract](physics.md) before selecting an accuracy label or
combining a decay response with a production calculation.

Version 0.26 retains the established numerical prescriptions. Its main changes
are package structure, scientific interfaces and redistributable examples.
