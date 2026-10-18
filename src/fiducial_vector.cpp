#include "leptonff/fiducial_vector.hpp"
#include "leptonff/evolution.hpp"
#include "leptonff/vector_current.hpp"
#include "leptonff/vector_decay.hpp"
#include "cut_breakpoints.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
Real pi() { return std::acos(-1.0L); }

Real cap_intersection(Real b1, Real b2, Real distance) {
  const Real a1 = std::acos(b1), a2 = std::acos(b2);
  if (distance >= a1+a2) return 0;
  if (distance <= std::abs(a1-a2)) return 2*pi()*(1-std::max(b1,b2));
  const Real sine = std::sin(distance);
  const Real one_minus_cosine = 2*std::pow(std::sin(distance/2),2);
  const Real cosine = 1-one_minus_cosine;
  const Real x = (b2-b1)+b1*one_minus_cosine;
  const Real determinant = std::max(Real{0},(1-b1)*(1+b1)*sine*sine-x*x);
  const Real root = std::sqrt(determinant);
  const Real A = std::atan2(root,x);
  const Real B = std::atan2(root,(b1-b2)+b2*one_minus_cosine);
  const Real C = std::atan2(root,cosine-b1*b2);
  return std::max(Real{0},2*(pi()-C-A*b1-B*b2));
}

Real band_acceptance_angle(Real b1, Real b2, Real angle) {
  if (b1 == 0 || b2 == 0) return 0;
  if (b1 == 1) return b2;
  if (b2 == 1) return b1;
  if (angle == 0) return std::min(b1,b2);
  return std::clamp(b1+b2-1+(cap_intersection(b1,b2,angle)
      +cap_intersection(b1,b2,pi()-angle))/(2*pi()),Real{0},std::min(b1,b2));
}

void validate_decay(Real mass, const FiducialCuts& cuts) {
  cuts.validate();
  MassWindow{cuts.mass_min,cuts.mass_max}.validate(mass);
}

Real angular_bound(Real momentum, const FiducialCuts& cuts) {
  if (momentum <= cuts.pt_min) return 0;
  const Real ratio = cuts.pt_min/momentum;
  return std::min(std::tanh(cuts.eta_max),std::sqrt((1-ratio)*(1+ratio)));
}

Real acceptance_unchecked(Real mass, Real muon, Real z, Real c,
                          const FiducialCuts& cuts,const SingleLeptonSelection& selection = {}) {
  const Real pair_mass = mass*std::sqrt(z);
  if (pair_mass <= cuts.mass_min || pair_mass >= cuts.mass_max) return 0;
  // c=1 is used as a collinear projection; massless endpoints are deliberately
  // excluded in the amplitude API, so construct these simple invariants here.
  const Real beta = std::sqrt(std::max(Real{0},1-4*std::pow(muon/mass,2)/z));
  const Real e1 = mass*((1+z)-(1-z)*beta*c)/4;
  const Real e2 = mass*((1+z)+(1-z)*beta*c)/4;
  const Real p1 = std::sqrt(std::max(Real{0},(e1-muon)*(e1+muon)));
  const Real p2 = std::sqrt(std::max(Real{0},(e2-muon)*(e2+muon)));
  const Real b1 = angular_bound(p1,cuts), b2 = angular_bound(p2,cuts);
  if (b1 == 0 || b2 == 0) return 0;
  // |p1 x p2| = |p1_transverse|*E_gamma. atan2 avoids losing the opening
  // angle when the pair is nearly back-to-back in the soft limit.
  const Real transverse = mass*std::sqrt(z)*beta*std::sqrt((1-c)*(1+c))/2;
  const Real cross = transverse*mass*(1-z)/2;
  const Real dot = e1*e2-(mass*mass*z-2*muon*muon)/2;
  const Real acute_angle = std::atan2(cross,std::abs(dot));
  if(selection.cut==SingleLeptonCut::none) return band_acceptance_angle(b1,b2,acute_angle);
  return (band_acceptance_angle(selection.restrict_bound(p1,b1),b2,acute_angle)
      +band_acceptance_angle(b1,selection.restrict_bound(p2,b2),acute_angle))/2;
}

std::vector<Real> outer_breaks(Real mass, const FiducialCuts& cuts,const SingleLeptonSelection& selection = {}) {
  auto result = one_leg_cut_breakpoints(mass,cuts,selection);
  result.push_back(std::max(Real{0},4*cuts.pt_min/mass-1));
  if(selection.cut==SingleLeptonCut::pt_above)
    for(Real t:{cuts.pt_min,selection.threshold}) result.push_back(4*t/mass-1);
  if(selection.cut==SingleLeptonCut::pt_above)
    result.push_back(2*(cuts.pt_min+selection.threshold)/mass-1);
  return result;
}

Real geometry_remainder(Real mass, Real muon, const FiducialCuts& cuts,
                         const HardDecayOptions& options,const SingleLeptonSelection& selection = {}) {
  if (cuts.pt_min == 0 && std::tanh(cuts.eta_max) == 1 && selection.cut==SingleLeptonCut::none) return 0;
  const Real threshold = 4*std::pow(muon/mass,2);
  const Real lower = std::max(threshold,std::pow(cuts.mass_min/mass,2));
  const Real upper = std::min(1-options.soft_cutoff,std::pow(cuts.mass_max/mass,2));
  if (lower >= upper) return 0;
  return integrate([&](Real z) {
    const Real reference = muon == 0 ? isotropic_acceptance(mass,z,1,cuts,selection)
        : acceptance_unchecked(mass,muon,z,1,cuts,selection);
    if (muon == 0) {
      const Real u_min = std::sqrt(options.collinear_cutoff);
      std::vector<Real> inner;
      if(selection.cut!=SingleLeptonCut::none) {
        std::vector<Real> energies{cuts.pt_min};
        if(std::tanh(cuts.eta_max)<1) energies.push_back(cuts.pt_min*std::cosh(cuts.eta_max));
        if(selection.cut==SingleLeptonCut::pt_above) {
          energies.push_back(selection.threshold);
          if(std::tanh(cuts.eta_max)<1) energies.push_back(selection.threshold*std::cosh(cuts.eta_max));
        } else if(std::tanh(selection.threshold)<1)
          energies.push_back(cuts.pt_min*std::cosh(selection.threshold));
        for(Real energy:energies) {
          const Real c=std::abs(((1+z)-4*energy/mass)/(1-z));
          if(c>0 && c<1) inner.push_back(std::sqrt(1-c));
        }
      }
      return integrate([&](Real u) {
        const Real c = 1-u*u;
        const Real difference = acceptance_unchecked(mass,0,z,c,cuts,selection)-reference;
        return 4*u*vector_scalar_real_density(mass,0,z,c)*difference;
      },u_min,1,options.integration,detail::resolved_cut_breakpoints(inner,u_min,1)).value;
    }
    const Real beta = std::sqrt(1-threshold/z);
    const Real edge = std::log((1+beta)/std::sqrt(threshold/z));
    return integrate([&](Real y) {
      const Real c = std::clamp(std::tanh(y)/beta,Real{0},Real{1});
      const Real jacobian = 1/(beta*std::pow(std::cosh(y),2));
      return 2*vector_scalar_real_density(mass,muon,z,c)*jacobian
          *(acceptance_unchecked(mass,muon,z,c,cuts,selection)-reference);
    },0,edge,options.integration).value;
  },lower,upper,options.integration,
      detail::resolved_cut_breakpoints(outer_breaks(mass,cuts,selection),lower,upper)).value;
}
} // namespace

Real rotational_band_acceptance(Real b1, Real b2, Real cosine) {
  if (!std::isfinite(b1) || !std::isfinite(b2) || b1 < 0 || b1 > 1 || b2 < 0 ||
      b2 > 1 || !std::isfinite(cosine) || std::abs(cosine) > 1)
    throw std::invalid_argument("Invalid rotational band geometry");
  return band_acceptance_angle(b1,b2,std::acos(std::abs(cosine)));
}

Real rotational_lepton_acceptance(Real p1,Real p2,Real cosine,
    const FiducialCuts& cuts,const SingleLeptonSelection& selection) {
  cuts.validate();selection.validate();
  if(!std::isfinite(p1) || !std::isfinite(p2) || p1<0 || p2<0)
    throw std::invalid_argument("Invalid observed lepton momentum");
  const Real b1=angular_bound(p1,cuts),b2=angular_bound(p2,cuts);
  if(selection.cut==SingleLeptonCut::none) return rotational_band_acceptance(b1,b2,cosine);
  return (rotational_band_acceptance(selection.restrict_bound(p1,b1),b2,cosine)
      +rotational_band_acceptance(b1,selection.restrict_bound(p2,b2),cosine))/2;
}

Real vector_radiative_acceptance(Real mass, Real muon, Real z, Real c,
                                 const FiducialCuts& cuts) {
  validate_decay(mass,cuts);
  if (!std::isfinite(muon) || muon < 0 || muon >= mass/2 || !std::isfinite(z) ||
      z <= 4*std::pow(muon/mass,2) || z >= 1 || !std::isfinite(c) || std::abs(c) > 1)
    throw std::invalid_argument("Invalid radiative acceptance point");
  return acceptance_unchecked(mass,muon,z,c,cuts);
}

Real vector_fragmented_acceptance(Real mass, Real z, Real c,
    Real f1, Real f2, const FiducialCuts& cuts,const SingleLeptonSelection& selection) {
  validate_decay(mass,cuts);
  selection.validate();
  if (!std::isfinite(z) || z <= 0 || z >= 1 || !std::isfinite(c) ||
      std::abs(c) > 1 || !std::isfinite(f1) || f1 < 0 || f1 > 1 ||
      !std::isfinite(f2) || f2 < 0 || f2 > 1)
    throw std::invalid_argument("Invalid fragmented vector point");
  const Real observed_mass = mass*std::sqrt(z*f1*f2);
  if (observed_mass <= cuts.mass_min || observed_mass >= cuts.mass_max) return 0;
  const Real e1 = mass*((1+z)-(1-z)*c)/4;
  const Real e2 = mass*((1+z)+(1-z)*c)/4;
  const Real b1 = angular_bound(e1*f1,cuts), b2 = angular_bound(e2*f2,cuts);
  if (b1 == 0 || b2 == 0) return 0;
  const Real cross = mass*mass*std::sqrt(z)*std::sqrt((1-c)*(1+c))*(1-z)/4;
  const Real dot = e1*e2-mass*mass*z/2;
  const Real angle=std::atan2(cross,std::abs(dot));
  if(selection.cut==SingleLeptonCut::none) return band_acceptance_angle(b1,b2,angle);
  return (band_acceptance_angle(selection.restrict_bound(e1*f1,b1),b2,angle)
      +band_acceptance_angle(b1,selection.restrict_bound(e2*f2,b2),angle))/2;
}

Real vector_scalar_real_density(Real mass, Real muon, Real z, Real c) {
  if (!std::isfinite(mass) || mass <= 0 || !std::isfinite(muon) || muon < 0 ||
      muon >= mass/2 || !std::isfinite(z) || z <= 4*std::pow(muon/mass,2) ||
      z >= 1 || !std::isfinite(c) || std::abs(c) > 1 || (muon == 0 && std::abs(c) == 1))
    throw std::invalid_argument("Invalid scalar real-emission point");
  const Real threshold = 4*std::pow(muon/mass,2);
  const Real beta2 = 1-threshold, pair_beta = std::sqrt(1-threshold/z);
  const Real v = 1-z, t = pair_beta*c;
  const Real denominator = (1-c)*(1+c)+threshold/z*c*c;
  // hep-ph/0505236v2 eq. (3), with its energy-fraction Jacobian. Keep the
  // finite-mass double pole; dropping it loses a finite collinear term.
  const Real A = -(3-beta2)*threshold/2;
  const Real B = (3-beta2)*(1+beta2)/v-2*(3-beta2)+2*v;
  return pair_beta*(2*A*(1+t*t)/(v*denominator*denominator)+B/denominator-v);
}

void HardDecayOptions::validate() const {
  if (!std::isfinite(soft_cutoff) || soft_cutoff <= 0 || soft_cutoff >= .1L ||
      !std::isfinite(collinear_cutoff) || collinear_cutoff <= 0 || collinear_cutoff >= .1L)
    throw std::invalid_argument("Hard remainder cutoffs must be in (0,0.1)");
  // Validate quadrature controls even for measurements with vanishing remainder.
  integrate([](Real) { return Real{0}; },0,0,integration);
}

FiducialVectorOrder1 fiducial_vector_order1(Real mass, const FiducialCuts& cuts,
    const QEDParameters& q, const HardDecayOptions& options,
    const IntegrationOptions& distribution_options,const SingleLeptonSelection& selection) {
  validate_decay(mass,cuts);
  q.validate();
  if (q.muon_mass >= mass/2) throw std::invalid_argument("Vector current below pair threshold");
  options.validate();
  selection.validate();
  FiducialVectorOrder1 result{};
  result.parent_mass = mass;
  result.cuts = cuts;
  result.parameters = q;
  result.born = isotropic_acceptance(mass,1,1,cuts,selection);
  result.projected_hard_coefficient = vector_mass_hard_coefficient(mass,q.mu_frag)
      .apply([&](Real z) { return isotropic_acceptance(mass,z,1,cuts,selection); },
             0,distribution_options,one_leg_cut_breakpoints(mass,cuts,selection));
  result.geometry_coefficient = geometry_remainder(mass,0,cuts,options,selection);
  result.hard_coefficient = result.projected_hard_coefficient+result.geometry_coefficient;
  result.ff_coefficient = decay_ff_coefficient_order1(mass,cuts,q,distribution_options,selection);
  result.full_coefficient = result.hard_coefficient+result.ff_coefficient;
  result.rate = result.born+q.expansion_parameter()*result.full_coefficient;
  return result;
}

MassiveFiducialVector massive_fiducial_vector_order1(Real mass,
    const FiducialCuts& cuts, const QEDParameters& q, const HardDecayOptions& options,
    const SingleLeptonSelection& selection) {
  validate_decay(mass,cuts);
  q.validate();
  options.validate();
  selection.validate();
  const Real muon = q.muon_mass;
  const Real born_ratio = massive_vector_born_ratio(mass,muon);
  const Real r4 = 4*std::pow(muon/mass,2);
  const Real momentum=mass*std::sqrt(1-r4)/2;
  const Real born_angular = selection.restrict_bound(momentum,angular_bound(momentum,cuts))
      *MassWindow{cuts.mass_min,cuts.mass_max}.born_measurement(mass);
  MassiveFiducialVector result{};
  result.born = born_ratio*born_angular;
  const Real projected = integrate([&](Real z) {
    return massive_vector_real_spectrum(mass,muon,z)
        *(acceptance_unchecked(mass,muon,z,1,cuts,selection)-born_angular);
  },r4,1-options.soft_cutoff,options.integration,outer_breaks(mass,cuts,selection)).value;
  result.coefficient = born_angular*massive_vector_inclusive_coefficient(mass,muon)
      +projected+geometry_remainder(mass,muon,cuts,options,selection);
  result.rate = result.born+q.expansion_parameter()*result.coefficient;
  return result;
}

FiducialVectorMatch fiducial_vector_ll_match(Real mass, const FiducialCuts& cuts,
    const LLBoundaryNonSinglet& evolution, const FiducialVectorOrder1& fixed,
    const IntegrationOptions& options) {
  const auto& q = evolution.parameters();
  const auto& p = fixed.parameters;
  if (mass != fixed.parent_mass || cuts.pt_min != fixed.cuts.pt_min ||
      cuts.eta_max != fixed.cuts.eta_max || cuts.mass_min != fixed.cuts.mass_min ||
      cuts.mass_max != fixed.cuts.mass_max || q.alpha != p.alpha ||
      q.muon_mass != p.muon_mass || q.mu_0 != p.mu_0 || q.mu_frag != p.mu_frag)
    throw std::invalid_argument("Hard coefficient and FF matching configurations differ");
  const auto ff = decay_ff_ll_match(mass,cuts,evolution,options);
  const Real addition = q.expansion_parameter()*fixed.hard_coefficient;
  return {ff.matched,addition,ff.matched+addition,fixed.rate};
}
} // namespace leptonff
