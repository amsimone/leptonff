# Photon-parent and matched pair-moment reference

Version 0.10.0, Release, AppleClang 21, 53-bit Real. Closed one-family MSbar
model, pure unpolarised vector at rest, massless vector Born normalisation.
M=91.1876 GeV, m=mu_0=0.1056583755 GeV. The explicit MSbar input at mu_0 is
1/137.035999084, without an on-shell conversion.

```sh
./build/leptonff_photon_pair configs/reference.conf
./build/leptonff_pair_moments configs/reference.conf
python3 python/check_photon_pair_reference.py build/leptonff_photon_pair \
  configs/reference.conf --output results/photon_pair_reference.json
python3 python/check_pair_moment_reference.py build/leptonff_pair_moments \
  configs/reference.conf --output results/pair_moment_reference.json
```

Saved results include source/executable/configuration hashes. The
[matching contract](../docs/photon_pair_matching.md) defines coefficients,
pair counting, coupling and limitations.

## Fixed-fraction fiducial photon-parent coefficients

pT>25 GeV, |eta|<2.4, 60<m_pair<120 GeV, mu_frag=M. These are coefficients
of alpha/(2*pi) for **one primary charge**, before photon-conversion FFs.
They are not rates.

| Photon fraction | Muon fraction | Projected | Finite geometry | Total |
| --- | --- | --- | --- | --- |
| 1 | 1 | -1.4430784883 | 0.3255797405 | -1.1174987478 |
| 0.9 | 0.95 | -1.1593164420 | 0.2572920830 | -0.9020243589 |
| 0.75 | 0.8 | -0.6918375608 | 0.0830683029 | -0.6087692579 |

Final sliver: t/(1-y)<2.5e-9. Changing from 1e-8 while tightening quadrature
shifts totals by at most 2.42e-8. Independent 96-node azimuth/energy checks
differ by at most 3.95e-7 in geometry and 2.21e-8 in projection; their 48-to-96
refinement is 2.32e-6. This supports roughly a few 1e-6 in the coefficient,
not every printed digit. The massive-current integral residual is 1.43e-13;
the analytic hard-moment residual is 7.57e-11. C++ tests also establish
hard+FF scale cancellation and decreasing mass-power differences.

## Inclusive hard-matched pair-mass moments

Weight each opposite-charge distinct-parent pair by `(m_pair²/M²)^(N-1)`.
No fiducial cuts. Two-loop running, mu_frag=M, hard alpha_MSbar(M):

| N | FF term | Muon hard term | Photon hard term | Matched moment |
| --- | --- | --- | --- | --- |
| 2 | 0.954824426476 | 0.011393736222 | -0.000051593558 | 0.966166569141 |
| 3 | 0.928571530405 | 0.017301874351 | -0.000016152950 | 0.945857251806 |
| 5 | 0.894759885481 | 0.025567616895 | -0.000005546257 | 0.920321956120 |

| N | mu_frag=M/2 | mu_frag=M | mu_frag=2M |
| --- | --- | --- | --- |
| 2 | 0.966193737994 | 0.966166569141 | 0.966115099477 |
| 3 | 0.945931486873 | 0.945857251806 | 0.945733634324 |
| 5 | 0.920504237737 | 0.920321956120 | 0.920044995011 |

Independent analytic-running/Magnus plus parent-tensor contraction agrees with
individual contributions/full results to 4.12e-13. Hard moments agree to
4.04e-10, expansion coefficients to 2.87e-8, hard coupling to 1.31e-18.
These are numerical checks, not physics uncertainties. The scale scan includes
beyond-NLL dependence from missing H2/d2 constants; it does not establish
per-mille accuracy for physical DY.

This joins hard and NLO FF ingredients for inclusive moments. The coupled
fiducial x-space hard convolution, physical spin/flavour content and NNLO host
remain separate deliverables. All 37 CTest checks pass in the reference build.
