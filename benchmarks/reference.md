# First fixed-order reference

Use the build/check commands in [README.md](../README.md) and
[reference.conf](../configs/reference.conf). The result is an isotropic, massless
parent-at-rest decay normalised to one, with the FF contribution only.

| Quantity | Initial reference |
| --- | ---: |
| Born fiducial acceptance | 0.836268494675 |
| FF coefficient multiplying α/(2π), both legs | −34.5761726595 |
| Born plus O(α) FF contribution | 0.796111391223 |

These are reference values, not a physical DY cross section. Small differences
in final digits can arise from platform precision and integration tolerances.
The output contains a tenfold tolerance-refinement comparison and per-parent
momentum residuals; the independent Python calculation checks the coefficient
to absolute 10⁻⁷ and records source, binary and configuration hashes. This check
tolerance corresponds to about 1.2 × 10⁻¹⁰ in the displayed rate.

The subtracted hard O(α) decay has not been added. Thus the numerical FF change
and fragmentation-scale dependence cannot yet be interpreted as a matched EW
correction or a perturbative uncertainty estimate.
