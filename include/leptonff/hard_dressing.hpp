#pragma once

#include "leptonff/fiducial_vector.hpp"
#include "leptonff/evolution.hpp"

namespace leptonff {
struct WeightedFraction { Real fraction, weight; };
// Symmetric finite muon-pair geometry at fixed fragmentation fractions.
// Shared by the scalar LL and coupled fiducial matchers.
Real vector_muon_pair_geometry(Real mass, const FiducialCuts& cuts,
    Real fraction1, Real fraction2, const HardDecayOptions& options = {},
    const SingleLeptonSelection& selection = {});
struct LLQuadrature {
  std::vector<WeightedFraction> nodes;
  Real number_residual = 0;
  Real first_moment_residual = 0;
  Real second_moment_residual = 0;
  Real third_moment_residual = 0;
  Real negative_weight_mass = 0;
};
// Retain the endpoint atom exactly. Panel the positive resolved measure in
// sqrt(1-z), using a two-node Gauss rule (moments through z^3) per panel. Preserve signed
// roundoff in a separate node; no clipping or renormalisation.
LLQuadrature compress_ll_quadrature(const LLNonSinglet& evolution,
                                    std::size_t panels);

struct HardDressingOptions {
  std::size_t ff_panels = 16;
  HardDecayOptions geometry;
  void validate() const;
};

struct HardProjection {
  std::vector<Real> initial_weights;
  std::vector<Real> evolved_weights;
  Real number_residual = 0;
  Real imaginary_weight_mass = 0;
  Real wrap_coefficient_bound = 0;
  Real omitted_tail_coefficient_bound = 0;
};
// H1 mass-spectrum distribution on the LL lattice, evolved on one leg.
// Its whole-plus endpoint is retained; the other leg is supplied at measurement.
HardProjection evolve_hard_projection(Real parent_mass, Real mu_frag,
                                      const LLNonSinglet& evolution);

struct HardDressedVectorMatch {
  Real ff_matched;
  Real undressed_hard_coefficient;
  Real projected_dressing_shift;
  Real geometry_dressing_shift;
  Real dressed_hard_coefficient;
  Real matched;
  Real additive_matched;
  HardProjection projection;
  LLQuadrature quadrature;
};

// LL on BOTH hard-event muons, with at most one finite insertion globally:
// UxU + a*(d1xU)xU + a*Ux(d1xU) + a*H1[UxU]. No H1*d1 term.
// Subtract the actual lattice FF overlap and restore the continuum O(alpha)
// coefficient, including the hard zero-evolution limit. Still lacks P1/running.
HardDressedVectorMatch fiducial_vector_hard_dressed_match(Real parent_mass,
    const FiducialCuts& cuts, const LLBoundaryNonSinglet& evolution,
    const FiducialVectorOrder1& fixed, const HardDressingOptions& options = {},
    const IntegrationOptions& distribution_options = {});
} // namespace leptonff
