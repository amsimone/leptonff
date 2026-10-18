# Canonical one-family benchmark

`order1.json` is a small, redistributable analytic reference. Inputs and the
absolute comparison tolerance are stored alongside its values. It does not
depend on MCFM, PHOTOS, a PDF set or excluded production archives.

Let `a=alpha/(2*pi)` and `L=log(mu_frag^2/muon_mass^2)`, with `mu_0=muon_mass`.
The full O(alpha) same-charge moments are:

- N=1: `1`.
- N=2: `1+a*(-4*L/3-17/9)`.
- N=3: `1+a*(-25*L/12-265/72)`.

For `P(z)=[(1+z^2)/(1-z)]_+`, an interior primitive is
`F(z)=-2*log(1-z)-z-z^2/2`. The interior bin is `F(.7)-F(.2)`;
the endpoint bin is `-F(.7)`. The convolution against a constant at `x=.3`
is `-log(x)+1/2+x+2*log(1-x)`.

The executable `build/examples/leptonff_benchmark` evaluates these quantities
through the numerical library. `tests/check_benchmark.py` compares its JSON to
the fixture; CTest registers this as `canonical_benchmark`.

These analytic values test definitions and endpoint actions. Full NLL evolution,
thresholds and matched W/Z calculations have separate invariant, independent
reference and regression checks; this fixture does not replace them.
