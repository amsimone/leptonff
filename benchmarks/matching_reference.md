# LL plus O(α) boundary, matched to the full O(α) FF contribution

```sh
./build/leptonff_match configs/reference.conf
python3 python/check_matching_reference.py build/leptonff_match configs/reference.conf \
  --output results/matching_reference.json
```

The configuration and isotropic parent-at-rest geometry match the earlier
references. Unlike the pure LL driver, this calculation evolves the finite
diagonal FF boundary with a single insertion across both muons. It then restores
the continuum O(α) FF coefficient by complete-overlap additive matching.

| Prediction | Dimensionless acceptance |
| --- | ---: |
| Full O(α) FF contribution on Born | 0.796111391223 |
| Matched FF-only, 262144 cells, epsilon=0.0002 | 0.796745916403 |
| Matched FF-only, 524288 cells, epsilon=0.0002 | 0.796745907742 |
| Matched FF-only, 524288 cells, epsilon=0.0001 | 0.796746821967 |

All domains have log span 32. The finer-grid shift is −8.7 × 10⁻⁹ and the next
cutoff shift is +9.1 × 10⁻⁷. The last row is not a zero-cutoff extrapolation.
The JSON records the LL rate, evolved-boundary coefficient, lattice overlap,
continuum fixed-order coefficient, additive correction and final matched rate.
The Python runner adds independent signed-convolution checks and source/binary/
configuration hashes.

This remains an FF-only reference. The NLO evolution kernel and hard O(α) decay
are absent, along with running/mixing and hadronic production. The difference
from fixed order is meaningful within this restricted prescription; it is not
a prediction of the complete NLL or physical NLO EW correction.
