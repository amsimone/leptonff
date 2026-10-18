#include "leptonff/nlo_qed.hpp"
#include "leptonff/vector_decay.hpp"

#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
// Compact QED expressions derived in docs/nlo_photon_kernels.md.
// Independent check: space-like QED (1606.02887v3) plus the complete
// analytic-continuation correction in hep-ph/9612250v1, eqs. (34)--(39).
Real lepton_to_photon(Real z) {
  const Real l = std::log(z), s = std::log1p(-z);
  const Real p = (1+(1-z)*(1-z))/z;
  return -.5L+4.5L*z+(-8+z/2)*l+2*z*s+(1-z/2)*l*l
      +p*(s*s-4*l*s-8*real_dilogarithm(z));
}

Real photon_to_lepton_direct(Real z) {
  const Real l = std::log(z), s = std::log1p(-z);
  const Real p = z*z+(1-z)*(1-z);
  const Real zeta2 = std::pow(std::acos(-1.0L),2)/6;
  const Real direct = -6+11.5L*z-10*z*z+(-2.5L+2*z+2*z*z)*l
      +(-3+2*z-2*z*z)*s+(-.5L+z-2*z*z)*l*l
      +p*(-2*zeta2+8*real_dilogarithm(z)+6*l*s-s*s);
  return direct;
}
Real photon_to_lepton_vacuum(Real z) {
  const Real l=std::log(z),s=std::log1p(-z),p=z*z+(1-z)*(1-z);
  return -4.0L/3+p*(-8.0L/9-4*(l+s)/3);
}
Real photon_to_lepton(Real z) {
  return photon_to_lepton_direct(z)+photon_to_lepton_vacuum(z);
}

Real photon_to_photon(Real z) {
  const Real l = std::log(z);
  return -4+12*z-164*z*z/9+92/(9*z)
      +(10+14*z+16*z*z/3+16/(3*z))*l+2*(1+z)*l*l;
}
} // namespace

Distribution nlo_photon_component(NLOPhotonComponent component) {
  switch(component) {
    case NLOPhotonComponent::fermion_to_photon:return {0,lepton_to_photon};
    case NLOPhotonComponent::photon_to_fermion_direct:return {0,photon_to_lepton_direct};
    case NLOPhotonComponent::photon_to_fermion_vacuum:return {0,photon_to_lepton_vacuum};
    case NLOPhotonComponent::photon_to_photon:return {-1,photon_to_photon,{},1};
  }
  throw std::invalid_argument("Unknown NLO photon component");
}
Distribution next_to_leading_kernel(Species parent, Species observed) {
  species_name(parent);
  species_name(observed);
  if (parent != Species::photon) {
    if (observed != Species::photon) return nlo_lepton_kernel(parent,observed);
    // Unlike the LO soft photon kernel, this coefficient is integrable at z=0.
    return {0,lepton_to_photon};
  }
  if (observed != Species::photon) return {0,photon_to_lepton}; // per charge
  // Abelian endpoint from the virtual two-loop photon term; no plus term.
  return {-1,photon_to_photon,{},1};
}
} // namespace leptonff
