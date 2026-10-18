#include "leptonff/fragmentation.hpp"

#include <cmath>
#include <stdexcept>

namespace leptonff {
const char* species_name(Species species) {
  switch (species) {
  case Species::muon_minus: return "mu-";
  case Species::muon_plus: return "mu+";
  case Species::photon: return "gamma";
  }
  throw std::invalid_argument("Unknown species");
}

void QEDParameters::validate() const {
  if (!std::isfinite(alpha) || alpha < 0 || alpha >= 1)
    throw std::invalid_argument("alpha must be finite and in [0,1)");
  for (const Real scale : {muon_mass, mu_0, mu_frag})
    if (!std::isfinite(scale) || scale <= 0)
      throw std::invalid_argument("QED masses/scales must be positive and finite");
}

Real QEDParameters::expansion_parameter() const {
  validate();
  return alpha / (2 * std::acos(-1.0L));
}

Distribution initial_condition_coefficient(Species parent, Species observed,
                                            const QEDParameters& parameters) {
  parameters.validate();
  species_name(parent);
  species_name(observed);
  const Real L = 2 * (std::log(parameters.mu_0) - std::log(parameters.muon_mass));
  if (parent != Species::photon) {
    if (observed == parent)
      return {0, {}, {[L](Real z) {
        return (1 + z*z) / (1-z) * (L - 2*std::log1p(-z) - 1);
      }}};
    if (observed == Species::photon)
      return {0, [L](Real z) {
        return (1 + (1-z)*(1-z)) / z * (L - 2*std::log(z) - 1);
      }, {}, 1}; // soft photon number moment diverges; momentum is finite
    return {}; // opposite-charge lepton first appears beyond O(alpha)
  }
  if (observed == Species::photon) return {-2 * L / 3};
  return {0, [L](Real z) { return (z*z + (1-z)*(1-z)) * L; }};
}

Distribution leading_kernel(Species parent, Species observed) {
  species_name(parent);
  species_name(observed);
  if (parent != Species::photon) {
    if (observed == parent)
      return {0, {}, {[](Real z) { return (1+z*z)/(1-z); }}};
    if (observed == Species::photon)
      return {0, [](Real z) { return (1+(1-z)*(1-z))/z; }, {}, 1};
    return {};
  }
  if (observed == Species::photon) return {-2.0L/3};
  return {0, [](Real z) { return z*z + (1-z)*(1-z); }};
}

Distribution fragmentation_coefficient_order1(Species parent, Species observed,
                                              const QEDParameters& parameters) {
  parameters.validate();
  const Real L = 2 * (std::log(parameters.mu_frag) - std::log(parameters.mu_0));
  return initial_condition_coefficient(parent, observed, parameters) +
         leading_kernel(parent, observed).scaled(L);
}

Distribution fragmentation_order1(Species parent, Species observed,
                                   const QEDParameters& parameters) {
  return Distribution(parent == observed ? 1 : 0) +
      fragmentation_coefficient_order1(parent, observed, parameters)
          .scaled(parameters.expansion_parameter());
}
} // namespace leptonff
