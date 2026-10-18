# Fixed-coupling photonic/non-singlet LL reference

Use [reference.conf](../configs/reference.conf). Both muons evolve independently
from a delta boundary at mu_0=m_mu, with fixed α and
`t = alpha/(2*pi) log(mu_frag²/mu_0²) = 0.0157033354839`.

```sh
./build/leptonff_ll configs/reference.conf
python3 python/check_ll_reference.py build/leptonff_ll configs/reference.conf \
  --output results/ll_reference.json
```

| Prediction | Dimensionless isotropic acceptance |
| --- | ---: |
| Born | 0.836268494675 |
| Consistently expanded LL O(α) | 0.798189960752 |
| LL, 262144 cells, epsilon=0.0002 | 0.798789753085 |
| LL, 524288 cells, epsilon=0.0002 | 0.798789799416 |
| LL, 524288 cells, epsilon=0.0001 | 0.798787837717 |

All lattices have log span 32. Grid refinement changes the rate by approximately
4.6 × 10⁻⁸; subsequent cutoff refinement changes it by −2.0 × 10⁻⁶. The last row
retains finite-cutoff effects and is not a zero-cutoff extrapolation. The finite
span wrap bound is about 3.6 × 10⁻¹³ per leg. Moment residuals and FFT roundoff
diagnostics are included in the JSON, as are source/executable/configuration
hashes when using the Python runner. Final digits vary with platform/tolerances.

The independent Poisson-series check uses a coarse lattice with separately
integrated rates, explicit virtual/no-emission weight, and direct convolutions;
its role is algorithm validation rather than a second fine-grid rate estimate.

This is a reference for the restricted LL model. It includes neither finite
O(α) FF boundaries nor hard O(α) decay, running α, photon conversion, PDFs or
QCD production. The older fixed-order FF benchmark includes finite FF terms;
differences between that benchmark and this LL result do not isolate higher
orders. No physical precision or scale-uncertainty claim follows from this table.
