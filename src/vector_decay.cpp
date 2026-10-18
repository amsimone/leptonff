#include "leptonff/vector_decay.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace leptonff {
namespace {
void positive_scale(Real scale) {
  if (!std::isfinite(scale) || scale <= 0)
    throw std::invalid_argument("Vector-decay masses/scales must be positive and finite");
}

struct MassiveKinematics {
  Real r, beta2, beta, t, logt, log_half, log_beta;
};

MassiveKinematics massive_kinematics(Real mass, Real muon) {
  positive_scale(mass);
  positive_scale(muon);
  if (muon >= mass/2)
    throw std::invalid_argument("Vector current must be above pair threshold");
  const Real r = std::pow(muon/mass,2);
  if (r == 0) throw std::invalid_argument("Muon mass ratio underflows");
  const Real beta2 = 1-4*r;
  const Real beta = std::sqrt(beta2);
  const Real t = 4*r/std::pow(1+beta,2); // stable (1-beta)/(1+beta)
  return {r,beta2,beta,t,std::log(t),std::log1p(-2*r/(1+beta)),
          std::log1p(-4*r)/2};
}
} // namespace

Distribution vector_mass_hard_coefficient(Real mass, Real mu_frag) {
  positive_scale(mass);
  positive_scale(mu_frag);
  const Real log_scale = 2*(std::log(mass)-std::log(mu_frag));
  return {1.5L,{}, {[log_scale](Real z) {
    return 2*(1+z*z)/(1-z)*(log_scale+std::log(z)+2*std::log1p(-z));
  }}};
}

Distribution vector_mass_leading_power_coefficient(Real mass, Real muon) {
  positive_scale(mass);
  positive_scale(muon);
  if (muon >= mass/2) throw std::invalid_argument("Vector current below pair threshold");
  const Real L = 2*(std::log(mass)-std::log(muon));
  return {1.5L,{}, {[L](Real z) {
    return 2*(1+z*z)/(1-z)*(L+std::log(z)-1);
  }}};
}

Distribution vector_mass_ff_coefficient(const QEDParameters& parameters) {
  return fragmentation_coefficient_order1(Species::muon_minus,Species::muon_minus,
                                           parameters).scaled(2);
}

Real massive_vector_born_ratio(Real mass, Real muon) {
  const auto k = massive_kinematics(mass,muon);
  return k.beta*(1+2*k.r);
}

Real massive_vector_real_spectrum(Real mass, Real muon, Real z) {
  const auto k = massive_kinematics(mass,muon);
  if (!std::isfinite(z) || z < 0 || z >= 1)
    throw std::invalid_argument("Massive real spectrum requires z in [0,1)");
  if (z <= 4*k.r) return 0;
  const Real ratio = 4*k.r/z;
  const Real pair_beta = std::sqrt(1-ratio);
  const Real photon_fraction = 1-z;
  // Bystritskiy et al., hep-ph/0505236v2, eq. (5). Its prefactor divided
  // by massless Born is exactly alpha/(2*pi). Avoid 1-pair_beta at small m.
  const Real logarithm = 2*std::log1p(pair_beta)-std::log(ratio);
  const Real B = (1+k.beta2)*(3-k.beta2)/photon_fraction
      -2*(3-k.beta2)+2*photon_fraction;
  return B*logarithm-2*((3-k.beta2)*z/photon_fraction+photon_fraction)*pair_beta;
}

Real massive_vector_real_integral(Real mass, Real muon, Real lower, Real upper,
                                  const IntegrationOptions& options) {
  const auto k = massive_kinematics(mass,muon);
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower < 0 ||
      upper < lower || upper >= 1)
    throw std::invalid_argument("Massive real integral requires 0<=lower<=upper<1");
  const Real start = std::max(lower,4*k.r);
  if (upper <= start) return 0;
  // Quadratic map resolves the square-root pair threshold. Direct bin
  // integration also avoids subtracting nearly equal cumulatives.
  const Real span = upper-start;
  return integrate([=](Real u) {
    return massive_vector_real_spectrum(mass,muon,start+span*u*u)*2*span*u;
  },0,1,options).value;
}

Real massive_vector_real_cumulative(Real mass, Real muon, Real z,
                                    const IntegrationOptions& options) {
  return massive_vector_real_integral(mass,muon,0,z,options);
}

Real massive_vector_inclusive_coefficient(Real mass, Real muon) {
  const auto k = massive_kinematics(mass,muon);
  const Real A = (1+k.beta2)/k.beta;
  const Real poly = 33+22*k.beta2-7*k.beta2*k.beta2;
  const Real log_coefficient = -poly/(8*k.beta*(3-k.beta2));
  // Combine cancelling log(m^2/M^2) terms before evaluating the high-energy limit.
  const Real eta = A*(4*real_dilogarithm(k.t)+2*real_dilogarithm(-k.t)
          -k.logt*(3*k.log_half-2*k.log_beta))
      +(3+log_coefficient)*std::log(k.r)-2*log_coefficient*k.log_half-4*k.log_beta
      +3*(5-3*k.beta2)/(4*(3-k.beta2));
  return 2*massive_vector_born_ratio(mass,muon)*eta;
}

void MassWindow::validate(Real mass) const {
  positive_scale(mass);
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower < 0 || upper <= lower)
    throw std::invalid_argument("Invalid vector mass window");
  if (upper == mass)
    throw std::invalid_argument("Upper mass cut exactly on the Born endpoint is infrared unsafe");
}

Real MassWindow::born_measurement(Real mass) const {
  validate(mass);
  return mass>lower && mass<upper ? 1 : 0;
}

Real MassWindow::apply(const Distribution& coefficient, Real mass,
                       const IntegrationOptions& options) const {
  validate(mass);
  const Real low = std::pow(lower/mass,2);
  const Real high = std::pow(upper/mass,2);
  return coefficient.apply([=](Real z) { return z>low && z<high ? 1.0L : 0.0L; },
                           0,options,{low,high});
}

VectorMassResult vector_mass_window(Real mass, const MassWindow& window,
                                    const QEDParameters& parameters,
                                    const IntegrationOptions& options) {
  parameters.validate();
  window.validate(mass);
  VectorMassResult result{};
  result.born = window.born_measurement(mass);
  result.hard_coefficient = window.apply(vector_mass_hard_coefficient(mass,parameters.mu_frag),mass,options);
  result.ff_coefficient = window.apply(vector_mass_ff_coefficient(parameters),mass,options);
  result.leading_power_coefficient = window.apply(
      vector_mass_leading_power_coefficient(mass,parameters.muon_mass),mass,options);
  result.matched_coefficient = result.hard_coefficient+result.ff_coefficient;
  const Real a = parameters.expansion_parameter();
  result.leading_power_rate = result.born+a*result.matched_coefficient;
  result.massive_born = result.born*massive_vector_born_ratio(mass,parameters.muon_mass);
  const Real low = std::pow(window.lower/mass,2);
  const Real high = std::pow(window.upper/mass,2);
  if (low < 1) {
    result.massive_coefficient = high>=1
        ? massive_vector_inclusive_coefficient(mass,parameters.muon_mass)
            -massive_vector_real_cumulative(mass,parameters.muon_mass,low,options)
        : massive_vector_real_integral(mass,parameters.muon_mass,low,high,options);
  }
  result.massive_rate = result.massive_born+a*result.massive_coefficient;
  return result;
}
} // namespace leptonff
