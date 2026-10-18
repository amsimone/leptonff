#pragma once
#include "leptonff/coupled_decay.hpp"

namespace leptonff {
struct CoupledMassBin {
  Real lower = 0, upper = 0; // GeV
  CoupledDecayMatch integral; // bin INTEGRAL, not divided by bin width
};

struct CoupledMassSpectrumResult {
  std::vector<CoupledMassBin> bins;
  CoupledDecayMatch total, sum;
  Real matched_closure = 0, fixed_order_closure = 0;
  Real coefficient1_closure = 0, coefficient2_closure = 0;
};

// Inclusive opposite-charge pair MASS spectrum in the same closed reference
// theory as CoupledFiducialDecay. Edges must strictly increase and exactly
// partition the parent fiducial mass window. Bin integrals are cumulative
// differences on ONE shared FF quadrature/support, with all hard channels and
// full overlap matching. Resolved total log-loss cells are integrated
// conservatively across bin edges; the zero-loss atom remains exact.
// No event-by-event covariance or statistical errors.
// Interior bins use (lower,upper]; the final upper edge retains the outer
// fiducial strict cut. In particular, an interior edge equal to M assigns the
// Born atom and its soft neighbourhood to the bin below M (IR safe).
class CoupledMassSpectrum {
public:
  CoupledMassSpectrum(Real mass, const FiducialCuts& cuts,
      const std::vector<Real>& edges, const QEDXSpaceGrid& grid,
      const CoupledDecayOptions& options = {});
  CoupledMassSpectrumResult match(const QEDXSpaceEvolution& evolution,
      const IntegrationOptions& fixed_order_options = {}) const;
  const std::vector<Real>& edges() const { return edges_; }
  const std::vector<CoupledFiducialDecay>& cumulative_operators() const { return cumulative_; }
private:
  std::vector<Real> edges_;
  // Above the kinematic endpoint the cumulative is exactly zero, so no hard
  // operator is constructed. Nonzero operators correspond to successive edges.
  std::vector<CoupledFiducialDecay> cumulative_;
};
} // namespace leptonff
