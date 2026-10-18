#pragma once

#include "leptonff/decay.hpp"

namespace leptonff {
class LLBoundaryNonSinglet;

// Average over a common uniform rotation of two fixed directions. Each cut
// requires |cos(theta_i)| < bound_i, with bounds in [0,1].
Real rotational_band_acceptance(Real bound1, Real bound2, Real opening_cosine);
// Mean one-lepton cumulative, with original pair cuts. Arguments are absolute
// three-momenta (equal to energies for massless particles).
Real rotational_lepton_acceptance(Real momentum1, Real momentum2, Real opening_cosine,
    const FiducialCuts& cuts, const SingleLeptonSelection& selection = {});
Real vector_radiative_acceptance(Real mass, Real muon_mass, Real z, Real cosine,
                                 const FiducialCuts& cuts);
// Massless real-event muons rescaled independently by fractions in [0,1].
// Directions are unchanged; the observed pair mass is M*sqrt(z*f1*f2).
Real vector_fragmented_acceptance(Real mass, Real z, Real cosine,
    Real fraction1, Real fraction2, const FiducialCuts& cuts,
    const SingleLeptonSelection& selection = {});
// Exact scalar real density, independently checked against Dirac traces.
// Supports m=0 away from the collinear endpoints.
Real vector_scalar_real_density(Real mass, Real muon_mass, Real z, Real cosine);

struct HardDecayOptions {
  IntegrationOptions integration{2e-6L,2e-6L,8192};
  Real soft_cutoff = 1e-6L;       // omit z>1-soft_cutoff in finite remainder
  Real collinear_cutoff = 1e-6L;  // omit |cosine|>1-collinear_cutoff, massless only
  void validate() const;
};

struct FiducialVectorOrder1 {
  Real parent_mass;
  FiducialCuts cuts;
  QEDParameters parameters;
  Real born;
  Real projected_hard_coefficient;
  Real geometry_coefficient;
  Real hard_coefficient;
  Real ff_coefficient;
  Real full_coefficient;
  Real rate;
};

// Unpolarised pure-vector parent at rest. Symmetric muon cuts only.
// Explicit soft/collinear cutoffs regulate the FINITE geometry remainder;
// neither replaces MSbar subtraction. Require refinement before interpretation.
FiducialVectorOrder1 fiducial_vector_order1(Real mass, const FiducialCuts& cuts,
    const QEDParameters& parameters, const HardDecayOptions& hard_options = {},
    const IntegrationOptions& distribution_options = {},
    const SingleLeptonSelection& selection = {});

struct MassiveFiducialVector {
  Real born;
  Real coefficient;
  Real rate;
};
// Exact massive real current + inclusive virtual normalisation. The finite
// soft remainder uses hard_options.soft_cutoff; no angular cutoff is needed.
MassiveFiducialVector massive_fiducial_vector_order1(Real mass,
    const FiducialCuts& cuts, const QEDParameters& parameters,
    const HardDecayOptions& hard_options = {},
    const SingleLeptonSelection& selection = {});

struct FiducialVectorMatch {
  Real ff_matched;
  Real hard_addition;
  Real matched;
  Real fixed_order;
};
// Additive NLO+LL in this restricted model: the full lattice FF overlap is
// already removed by decay_ff_ll_match. H1 is inserted once, WITHOUT LL dressing.
// Does not claim NLL or the hard-generated subleading logarithmic towers.
FiducialVectorMatch fiducial_vector_ll_match(Real mass, const FiducialCuts& cuts,
    const LLBoundaryNonSinglet& evolution, const FiducialVectorOrder1& fixed,
    const IntegrationOptions& options = {});
} // namespace leptonff
