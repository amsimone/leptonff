# Generator-independent production/decay interface

`production_decay.hpp` is the common observable-level interface. It accepts
explicit configuration objects; it does not read environment variables or
depend on production PDFs, a shower or a host's histogram machinery.

`DecayCalculator` owns immutable evolution and matching tables. Its
`cumulative` and `evaluate` methods handle Born-decay points; `response`
returns decay tensors that a production density can contract. The MCFM
physical-spectrum callback now uses this common evaluator. Its historical
rate/histogram ABI is preserved as a compatibility adapter, not the public
library API. That adapter still has 17 online bins; arbitrary finer bins
are available through the common API and offline replay.

## Weighted Born-decay points

`BornDecayPoint` contains first and second massless four-vectors and a
64-bit radiation seed. First means the charged W lepton, or negative Z
lepton; second means the neutrino or positive Z lepton. Native library and
C-interface momenta are `(energy,px,py,pz)` in GeV. MCFM's `(px,py,pz,E)`
layout is translated only in its adapter. Tiny host mapping residuals may
be repaired by the explicit `map_massless_parent` helper with a stated
tolerance; generic physics functions do not silently repair arbitrary input.

The returned functional has Born, order-alpha, matched, second-order and
higher-order channels. Production weights are external. The muon/neutrino
interpretation, unitary branching-ratio convention, perturbative boundaries,
photonic conditioning and inclusive-pair definition are physics choices.
These are not a generic shower algorithm for arbitrary hard processes.

Finite W/Z matching reconstructs spin from the original scalar Born decay
through the dual angular projector. This requires full Born angular
coverage at fixed production kinematics, without original-lepton-dependent
cuts/scales. Original real and mapped subtraction terms must use the SAME
radiation stream. A weighted LHE file without these guarantees is not
automatically a valid input. QED FSR must be absent or explicitly removed.

## Direct production tensors

`SpinDecayResponse` contains five Hermitian 3x3 tensors per rate/bin.
Contract them with a signed production density using `contract`. The frame
is Cartesian axes obtained by a rotation-free boost from laboratory axes.
No division by a Born rate or positivity requirement is introduced.

`drell_yan.hpp` converts UNNORMALISED angular coefficients `(U,H0,...,H7)`
into a production tensor, with the complete convention stated in its header.
`Hi=U*Ai` must be supplied as weights, not evaluated by dividing an unstable
or zero matching contribution. W odd channels and both charges are tested
against the explicit angular polynomial. For Z only the parity-even part
is retained, valid for the charge-symmetric pair observables implemented here.
Forward-backward or charge-asymmetric Z measurements require an extension.

`collins_soper_rotation` defines an explicit beam-bisector frame, including
its x-axis sign and the qT=0 choice. A host must verify that convention and
its recoil prescription; the existence of this helper does NOT prove that
DYTurbo's native axes already agree with it.

The response uses Gauss-Legendre/trapezoidal Born-angle quadrature plus
separate finite-photon sampling. Both need refinement. Tests check the
inclusive Born tensor against the identity and the finite response against
an independent weighted angular-tomography integral. Returned tensors can
be retained and reused with new production weight variations at the SAME
phase-space node and decay configuration. A validated interpolating
response grid and its quadrature-error accounting remain future work.

## Production records and covariance

`production_records.hpp` provides a streaming, versioned text protocol:

```text
LEPTONFF_PRODUCTION 1
META process "producer" "physics_id" "run_id" "fb" random_seed 1 0 nweights
WEIGHTS "central" ...
GROUP id nterms
POINT seed first_E first_px first_py first_pz second_E ... weight ...
...
BLOCK part channel iteration trials stored_groups
...
END
```

Numbers are written with `long double` round-trip precision. Each group is
buffered until the generator accepts/commits it; rejected groups are dropped.
All weights and kinematics of a real/subtraction group are preserved.
Zero/unrepresented trials remain in the BLOCK denominator. END is mandatory;
interrupted archives, duplicate blocks, inconsistent seeds and invalid
coverage/QED metadata are rejected. Metadata/file hashes in the completed
run manifest are checked by the recommended replay front end.

For a group vector x, the block stores sums S and products P only AFTER
summing its signed terms. Host weights already include integration
normalisation, so the block mean is S, not S/ncall. Its covariance is
`(ncall*P-S*S^T)/(ncall-1)`. Equal iteration means and covariances are then
combined within each part/channel; different contributions are summed.
Each resumed iteration uses its own call count. Cross-model/channel
products are kept and cannot be reconstructed from marginal errors.

This estimator is for Monte Carlo integration groups. Deterministic
quadrature nodes are not independent random events and must not be fed
into it to invent Monte Carlo errors. A quadrature backend should contract
decay responses and carry its integration-error information separately.

MCFM production-only export bypasses fragmentation/PHOTOS entirely.
Its physics ID includes production settings, PDFs, host executable and
branching inputs while excluding output paths, integration seed/budget and
technical subtraction parameters. The separate run ID and random seed
identify the sampling. Archives with the SAME random seed are not
independent just because their output paths differ. Automatic pooling of
arbitrary archives is not supplied; independence must be established.
