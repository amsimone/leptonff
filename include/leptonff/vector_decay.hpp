#pragma once

#include "leptonff/fragmentation.hpp"
#include "leptonff/special_functions.hpp"

namespace leptonff {
// Angular-inclusive unpolarised VECTOR current only, z=m_mumu^2/M^2.
// Coefficients multiply a=alpha/(2*pi), normalised to massless vector Born.
// No pT/eta cuts or axial-current finite-mass terms are represented here.
Distribution vector_mass_hard_coefficient(Real parent_mass, Real mu_frag);
Distribution vector_mass_leading_power_coefficient(Real parent_mass, Real muon_mass);
Distribution vector_mass_ff_coefficient(const QEDParameters& parameters);


Real massive_vector_born_ratio(Real parent_mass, Real muon_mass);
// Exact massive real spectrum, coefficient dC/dz; zero below pair threshold.
Real massive_vector_real_spectrum(Real parent_mass, Real muon_mass, Real z);
Real massive_vector_real_integral(Real parent_mass, Real muon_mass, Real lower,
                                  Real upper, const IntegrationOptions& options = {});
// Exact massive hard-photon cumulative: coefficient on z in [4m^2/M^2,z_max].
// Upper endpoint z_max must be <1; virtual/soft terms are not part of this.
Real massive_vector_real_cumulative(Real parent_mass, Real muon_mass, Real z_max,
                                    const IntegrationOptions& options = {});
// Exact massive inclusive O(a) coefficient, including virtual+real radiation.
Real massive_vector_inclusive_coefficient(Real parent_mass, Real muon_mass);

struct MassWindow {
  Real lower = 60; // GeV, strict mass_min < m_mumu < mass_max
  Real upper = 120;
  void validate(Real parent_mass) const;
  Real born_measurement(Real parent_mass) const;
  Real apply(const Distribution& coefficient, Real parent_mass,
             const IntegrationOptions& options = {}) const;
};

struct VectorMassResult {
  Real born;
  Real hard_coefficient;
  Real ff_coefficient;
  Real leading_power_coefficient;
  Real matched_coefficient;
  Real leading_power_rate;
  Real massive_born;
  Real massive_coefficient;
  Real massive_rate;
};

// First-order hard+FF matching and an independent massive analytic reference.
VectorMassResult vector_mass_window(Real parent_mass, const MassWindow& window,
                                    const QEDParameters& parameters,
                                    const IntegrationOptions& options = {});
} // namespace leptonff
