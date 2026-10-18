#pragma once

#include "leptonff/fragmentation.hpp"

#include <utility>

namespace leptonff {
class LLNonSinglet;
class LLBoundaryNonSinglet;
struct FourVector {
  Real energy, px, py, pz;
  FourVector scaled(Real z) const;
  Real pt() const;
  Real eta() const;
};

Real invariant_mass(const FourVector& first, const FourVector& second);
std::pair<FourVector, FourVector> born_decay(Real mass, Real cos_theta, Real phi = 0);

struct FiducialCuts {
  Real pt_min = 25;       // GeV; both muons
  Real eta_max = 2.4L;
  Real mass_min = 60;     // GeV
  Real mass_max = 120;    // GeV
  void validate() const;
  bool accepts(const FourVector& first, const FourVector& second) const;
};

enum class SingleLeptonCut { none, pt_above, absolute_eta_below };
// Pair-weighted mean of applying a cumulative cut to either observed lepton:
// (I_1+I_2)/2, with BOTH original fiducial lepton cuts still imposed.
// A histogram counting both muons multiplies this measurement by two.
struct SingleLeptonSelection {
  SingleLeptonCut cut = SingleLeptonCut::none;
  Real threshold = 0;
  void validate() const;
  Real restrict_bound(Real momentum, Real fiducial_bound) const;
};

// Normalised isotropic, massless two-body decay of a parent at rest. Angular
// integration is analytic. This is a fixture, not the physical DY spin density.
Real isotropic_acceptance(Real mass, Real z1, Real z2, const FiducialCuts& cuts,
    const SingleLeptonSelection& selection = {});
std::vector<Real> one_leg_cut_breakpoints(Real mass, const FiducialCuts& cuts,
    const SingleLeptonSelection& selection = {});

// Coefficient of a in Born x [I x I + a(D1 x I + I x D1)]. Hard O(alpha)
// decay is absent: this is the FF contribution only, not matched NLO EW.
Real decay_ff_coefficient_order1(Real mass, const FiducialCuts& cuts,
                                 const QEDParameters& parameters,
                                 const IntegrationOptions& options = {},
                                 const SingleLeptonSelection& selection = {});

// Identical independent LL FFs on both legs. Prefix sums account for the
// coupled mass window without an O(bins^2) tensor loop.
Real decay_ll_acceptance(Real mass, const FiducialCuts& cuts,
                         const LLNonSinglet& evolution);

// Bilinear lattice measurement for two signed weight arrays on the same grid.
// If average_mass_cells, integrate each resolved total log-loss shell over its
// centred cell of width spacing. The zero-loss Born atom stays exact. Angular
// cuts are still evaluated at the grid fractions. This is a numerical bin
// integration prescription with a vanishing cell width, not detector smearing.
Real decay_lattice_acceptance(Real mass, const FiducialCuts& cuts, Real spacing,
                              const std::vector<Real>& first,
                              const std::vector<Real>& second,
                              bool average_mass_cells = false,
                              const SingleLeptonSelection& selection = {});

struct FFDecayMatch {
  Real born;
  Real ll;
  Real evolved_boundary_coefficient;
  Real resummed_with_boundary;
  Real lattice_overlap_coefficient;
  Real fixed_order_coefficient;
  Real additive_correction;
  Real matched;
};

// Matches ONLY the FF contribution through O(alpha), including the finite
// boundary and the lattice's exact first-order overlap. No hard decay yet.
FFDecayMatch decay_ff_ll_match(Real mass, const FiducialCuts& cuts,
                              const LLBoundaryNonSinglet& evolution,
                              const IntegrationOptions& options = {});
} // namespace leptonff
