#include "leptonff/nlo_lepton.hpp"
#include "leptonff/special_functions.hpp"

#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
constexpr Real zeta3 = 1.20205690315959428539973816151145L;
Real zeta2() { return std::pow(std::acos(-1.0L), 2)/6; }

// Arbuzov--Melnikov hep-ph/0205172v1, eqs. (15)--(18), with delta and
// -20/9 [1/(1-z)]_+ kept explicitly. No added whole-plus prescription.
Real photonic(Real z) {
  const Real l = std::log(z);
  return (1+z*z)/(1-z)*(2*l*std::log1p(-z)-2*l*l
      -2*real_dilogarithm(1-z)) + (1+z)*l*l/2 + 2*z*l-3*z+2;
}
Real nonsinglet(Real z) {
  return -2.0L/3*(1+z*z)*std::log(z)/(1-z)-2.0L/9+22*z/9;
}
Real singlet(Real z) {
  const Real l = std::log(z);
  return (1+z)*l*l+(-5-9*z-8*z*z/3)*l-8-20/(9*z)+4*z+56*z*z/9;
}
Real same_interference(Real z) {
  return (1+z*z)/(1-z)*(2*real_dilogarithm(1-z)+1.5L*std::log(z))
      -3.5L*(1+z)*std::log(z)-7+8*z;
}
// Arbuzov hep-ph/0206036v2, eq. (5.15). This is NOT half the same-charge
// function pointwise; only their first moments have the factor of two.
Real opposite_interference(Real z) {
  const Real l = std::log(z);
  return 2*(1+z*z)/(1+z)*(-2*real_dilogarithm(-z)
      -2*l*std::log1p(z)+l*l/2-zeta2())+2*(1+z)*l+4*(1-z);
}
// Algebraic sum photonic + same-charge interference cancels Li2(1-z).
// The non-singlet density agrees with Blumlein--Kawamura hep-ph/0211191v1
// eq. (30), when its singular p(z) term is split into regular/plus/delta.
Real diagonal_without_singlet(Real z) {
  const Real l = std::log(z);
  return (1+z*z)/(1-z)*(2*l*std::log1p(-z)-2*l*l+1.5L*l)
      +(1+z)*l*l/2-(7+3*z)*l/2-5*(1-z)+nonsinglet(z);
}
Real diagonal_delta() { return 3.0L/8-3*zeta2()+6*zeta3-4*zeta2()/3-1.0L/6; }
Real pair_plus(Real z) { return -20/(9*(1-z)); }
} // namespace

const char* nlo_lepton_component_name(NLOLeptonComponent component) {
  switch (component) {
  case NLOLeptonComponent::photonic: return "photonic";
  case NLOLeptonComponent::nonsinglet_pair: return "nonsinglet_pair";
  case NLOLeptonComponent::singlet_pair: return "singlet_pair";
  case NLOLeptonComponent::same_charge_interference: return "same_charge_interference";
  case NLOLeptonComponent::opposite_charge_interference: return "opposite_charge_interference";
  }
  throw std::invalid_argument("Unknown NLO lepton component");
}

Distribution nlo_lepton_component(NLOLeptonComponent component) {
  nlo_lepton_component_name(component);
  switch (component) {
  case NLOLeptonComponent::photonic:
    return {3.0L/8-3*zeta2()+6*zeta3, photonic};
  case NLOLeptonComponent::nonsinglet_pair:
    return {-4*zeta2()/3-1.0L/6, nonsinglet, {pair_plus}};
  case NLOLeptonComponent::singlet_pair:
    return {0, singlet, {}, 1}; // 1/z: number moment diverges
  case NLOLeptonComponent::same_charge_interference:
    return {0, same_interference};
  case NLOLeptonComponent::opposite_charge_interference:
    return {0, opposite_interference};
  }
  throw std::invalid_argument("Unknown NLO lepton component");
}

Distribution nlo_lepton_kernel(Species parent, Species observed) {
  species_name(parent);
  species_name(observed);
  if (parent == Species::photon || observed == Species::photon)
    throw std::invalid_argument("NLO photon channels have not been implemented");
  if (parent == observed)
    return {diagonal_delta(), [](Real z) {
      return diagonal_without_singlet(z)+singlet(z);
    }, {pair_plus}, 1};
  return {0, [](Real z) { return singlet(z)+opposite_interference(z); }, {}, 1};
}

Distribution nlo_valence_kernel() {
  return {diagonal_delta(), [](Real z) {
    return diagonal_without_singlet(z)-opposite_interference(z);
  }, {pair_plus}};
}
} // namespace leptonff
