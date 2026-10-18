#pragma once
#include "leptonff/moment_evolution.hpp"
#include "leptonff/decay.hpp"

namespace leptonff {
using LatticeMatrix = std::array<std::array<std::vector<Real>,3>,3>;

struct QEDXSpaceGrid {
  std::size_t bins = 4096; // power of two, u_j=j*log_span/bins
  Real log_span = 40;
  Real wrap_exponent = .75L; // exponential bound, strictly between 0 and 1
  Real max_wrap_bound = 1e-10L;
  Real max_tail_bound = 1e-10L;
  void validate() const;
  Real spacing() const { return log_span/bins; }
};

struct QEDXSpaceDiagnostics {
  Real weighted_wrap_bound = 0, weighted_tail_bound = 0;
  Real imaginary_mass = 0, negative_weight_mass = 0;
  Real maximum_frequency_step_shift = 0;
  Real weighted_step_shift_bound = 0;
  std::size_t maximum_steps = 0;
};

// Every lattice array represents z*D(z) dz, NOT D(z) dz. This energy tilt
// makes the full singlet/photon matrix integrable at small z. No clipping.
// Piecewise-linear interpolation in z resolves the endpoint plus/delta action
// and preserves the N=2 and N=3 kernel moments before tail omission.
class QEDXSpaceEvolution {
public:
  explicit QEDXSpaceEvolution(const MomentEvolutionParameters& parameters = {},
      const QEDXSpaceGrid& grid = {}, const MomentEvolutionOptions& options = {});
  const MomentEvolutionParameters& parameters() const { return parameters_; }
  const QEDXSpaceGrid& grid() const { return grid_; }
  const QEDXSpaceDiagnostics& diagnostics() const { return diagnostics_; }
  const LatticeMatrix& weights() const { return weights_; }
  const LatticeMatrix& propagator() const { return propagator_; }
  const LatticeMatrix& kernel0() const { return kernel0_; }
  const LatticeMatrix& kernel1() const { return kernel1_; }
  const LatticeMatrix& boundary1() const { return boundary1_; }
  const LatticeMatrix& coefficient1() const { return coefficient1_; }
  const LatticeMatrix& coefficient2() const { return coefficient2_; }
  const std::vector<std::array<std::vector<Real>,3>>& flavour_weights() const { return flavour_weights_; }
  std::vector<std::array<Real,3>> flavour_moment(Real N) const;
  Real z(std::size_t index) const;
  MomentMatrix moment(Real N) const; // N>=2: bounded on the energy-weighted grid
  // Truncated action on z>=minimum_z>0. Undo the tilt only on this support;
  // full inclusive number integrals need not exist. Lower-cut refinement is
  // a measurement/grid error, distinct from periodic-wrap/tail bounds.
  Real apply(Species parent, Species observed, const Function& measurement,
             Real minimum_z) const;
private:
  MomentEvolutionParameters parameters_;
  QEDXSpaceGrid grid_;
  QEDXSpaceDiagnostics diagnostics_;
  LatticeMatrix weights_, propagator_, kernel0_, kernel1_, boundary1_, coefficient1_, coefficient2_;
  std::vector<std::array<std::vector<Real>,3>> flavour_weights_;
};

struct XSpaceFFDecayMatch {
  Real born = 0, resummed = 0, charge_flipped_pairs = 0;
  Real lattice_coefficient1 = 0, continuum_coefficient1 = 0;
  Real additive_correction = 0, matched = 0;
};
// Born two-direction, inclusive identified opposite-charge pair density:
// D_same x D_same + D_opposite x D_opposite, one observed lepton per direction.
// Exactly O(alpha) FF matched, with no hard H1 and no exclusive extra-pair veto.
// This is FF-layer validation, NOT a full fiducial NLL cross section.
XSpaceFFDecayMatch decay_xspace_ff_match(Real mass,const FiducialCuts& cuts,
    const QEDXSpaceEvolution& evolution,const IntegrationOptions& options = {});
} // namespace leptonff
