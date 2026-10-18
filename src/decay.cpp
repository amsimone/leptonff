#include "leptonff/decay.hpp"
#include "leptonff/evolution.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace leptonff {
namespace {
void require_mass(Real mass) {
  if (!std::isfinite(mass) || mass <= 0)
    throw std::invalid_argument("Decay mass must be positive and finite");
}
void require_fraction(Real z) {
  if (!std::isfinite(z) || z < 0 || z > 1)
    throw std::invalid_argument("Collinear fraction must be in [0,1]");
}
} // namespace

FourVector FourVector::scaled(Real z) const {
  require_fraction(z);
  return {energy*z, px*z, py*z, pz*z};
}
Real FourVector::pt() const { return std::hypot(px, py); }
Real FourVector::eta() const {
  if (pt() == 0) {
    if (pz == 0) throw std::domain_error("Zero momentum has no pseudorapidity");
    return std::copysign(std::numeric_limits<Real>::infinity(), pz);
  }
  return std::asinh(pz / pt());
}

Real invariant_mass(const FourVector& first, const FourVector& second) {
  const Real E = first.energy + second.energy;
  const Real x = first.px + second.px;
  const Real y = first.py + second.py;
  const Real z = first.pz + second.pz;
  const Real square = E*E - x*x - y*y - z*z;
  if (square < -64*std::numeric_limits<Real>::epsilon()*E*E)
    throw std::domain_error("Spacelike pair in invariant_mass");
  return std::sqrt(std::max(Real{0}, square));
}

std::pair<FourVector, FourVector> born_decay(Real mass, Real cos_theta, Real phi) {
  require_mass(mass);
  if (!std::isfinite(cos_theta) || std::abs(cos_theta) > 1 || !std::isfinite(phi))
    throw std::invalid_argument("Invalid decay angles");
  const Real E = mass / 2;
  const Real pt = E * std::sqrt((1-cos_theta)*(1+cos_theta));
  return {{E, pt*std::cos(phi), pt*std::sin(phi), E*cos_theta},
          {E, -pt*std::cos(phi), -pt*std::sin(phi), -E*cos_theta}};
}

void FiducialCuts::validate() const {
  for (Real x : {pt_min, eta_max, mass_min, mass_max})
    if (!std::isfinite(x)) throw std::invalid_argument("Non-finite fiducial cut");
  if (pt_min < 0 || eta_max <= 0 || mass_min < 0 || mass_max <= mass_min)
    throw std::invalid_argument("Invalid fiducial cut range");
}

bool FiducialCuts::accepts(const FourVector& first, const FourVector& second) const {
  validate();
  if (first.pt() <= pt_min || second.pt() <= pt_min) return false;
  const Real m = invariant_mass(first, second);
  return std::abs(first.eta()) < eta_max && std::abs(second.eta()) < eta_max &&
         m > mass_min && m < mass_max;
}

void SingleLeptonSelection::validate() const {
  if(!std::isfinite(threshold) || threshold<0 ||
      (cut!=SingleLeptonCut::none && cut!=SingleLeptonCut::pt_above &&
       cut!=SingleLeptonCut::absolute_eta_below))
    throw std::invalid_argument("Invalid single-lepton cumulative selection");
}
Real SingleLeptonSelection::restrict_bound(Real momentum,Real base) const {
  if(cut==SingleLeptonCut::none) return base;
  if(cut==SingleLeptonCut::absolute_eta_below) return std::min(base,std::tanh(threshold));
  if(momentum<=threshold) return 0;
  const Real r=threshold/momentum;
  return std::min(base,std::sqrt((1-r)*(1+r)));
}

Real isotropic_acceptance(Real mass, Real z1, Real z2, const FiducialCuts& cuts,
    const SingleLeptonSelection& selection) {
  require_mass(mass);
  require_fraction(z1);
  require_fraction(z2);
  cuts.validate();
  selection.validate();
  if (z1 == 0 || z2 == 0) return 0;
  const Real observed_mass = mass * std::sqrt(z1*z2);
  if (observed_mass <= cuts.mass_min || observed_mass >= cuts.mass_max) return 0;
  const Real ratio = 2*cuts.pt_min/(mass*std::min(z1, z2));
  if (ratio >= 1) return 0;
  const Real base=std::min(std::tanh(cuts.eta_max), std::sqrt((1-ratio)*(1+ratio)));
  return (selection.restrict_bound(mass*z1/2,base)+selection.restrict_bound(mass*z2/2,base))/2;
}

std::vector<Real> one_leg_cut_breakpoints(Real mass, const FiducialCuts& cuts,
    const SingleLeptonSelection& selection) {
  require_mass(mass);
  cuts.validate();
  selection.validate();
  // Includes both discontinuities (mass) and the pt/eta branch transition.
  std::vector<Real> points{2*cuts.pt_min/mass,
      std::pow(cuts.mass_min/mass, 2), std::pow(cuts.mass_max/mass, 2)};
  // Once tanh(eta_max) rounds to one there is no resolved eta branch; avoid
  // cosh overflow in an otherwise legitimate inclusive diagnostic.
  if (std::tanh(cuts.eta_max) < 1 && cuts.pt_min > 0)
    points.push_back(2*cuts.pt_min*std::cosh(cuts.eta_max)/mass);
  if(selection.cut==SingleLeptonCut::pt_above && selection.threshold>0) {
    points.push_back(2*selection.threshold/mass);
    points.push_back(cuts.pt_min/selection.threshold);
    if(std::tanh(cuts.eta_max)<1)
      points.push_back(2*selection.threshold*std::cosh(cuts.eta_max)/mass);
  }
  if(selection.cut==SingleLeptonCut::absolute_eta_below && std::tanh(selection.threshold)<1)
    points.push_back(2*cuts.pt_min*std::cosh(selection.threshold)/mass);
  return points;
}

Real decay_ff_coefficient_order1(Real mass, const FiducialCuts& cuts,
                                 const QEDParameters& parameters,
                                 const IntegrationOptions& options,
                                 const SingleLeptonSelection& selection) {
  const auto first = fragmentation_coefficient_order1(
      Species::muon_plus, Species::muon_plus, parameters);
  const auto second = fragmentation_coefficient_order1(
      Species::muon_minus, Species::muon_minus, parameters);
  const auto breaks = one_leg_cut_breakpoints(mass, cuts,selection);
  const auto measurement = [&](Real z1, Real z2) {
    return isotropic_acceptance(mass, z1, z2, cuts,selection);
  };
  // Do not multiply two truncated FFs: that would add an uncontrolled a^2 term.
  return first.apply([&](Real z) { return measurement(z, 1); }, 0, options, breaks) +
         second.apply([&](Real z) { return measurement(1, z); }, 0, options, breaks);
}

Real decay_ll_acceptance(Real mass, const FiducialCuts& cuts,
                         const LLNonSinglet& evolution) {
  return decay_lattice_acceptance(mass,cuts,evolution.grid().spacing(),
                                  evolution.weights(),evolution.weights());
}

Real decay_lattice_acceptance(Real mass, const FiducialCuts& cuts, Real h,
                              const std::vector<Real>& first_weights,
                              const std::vector<Real>& second_weights,
                              bool average_mass_cells,const SingleLeptonSelection& selection) {
  require_mass(mass);
  cuts.validate();
  selection.validate();
  if (first_weights.empty() || first_weights.size() != second_weights.size() ||
      !std::isfinite(h) || h <= 0)
    throw std::invalid_argument("Incompatible decay lattice weights/spacing");
  for (std::size_t j = 0; j < first_weights.size(); ++j)
    if (!std::isfinite(first_weights[j]) || !std::isfinite(second_weights[j]))
      throw std::invalid_argument("Non-finite decay lattice weight");
  const Real mass_lower = 2*(std::log(mass)-std::log(cuts.mass_max));
  const Real mass_upper = cuts.mass_min == 0 ? std::numeric_limits<Real>::infinity()
      : 2*(std::log(mass)-std::log(cuts.mass_min));
  const Real eta_bound = selection.cut==SingleLeptonCut::absolute_eta_below
      ? std::min(std::tanh(cuts.eta_max),std::tanh(selection.threshold)) : std::tanh(cuts.eta_max);
  if(selection.cut==SingleLeptonCut::pt_above && selection.threshold>cuts.pt_min) {
    std::vector<Real> base(first_weights.size(),0),tight(base);
    for(std::size_t j=0;j<base.size();++j) {
      const Real momentum=mass*std::exp(-h*j)/2;
      if(momentum<=cuts.pt_min) continue;
      const Real r=cuts.pt_min/momentum;
      base[j]=std::min(eta_bound,std::sqrt((1-r)*(1+r)));
      tight[j]=selection.restrict_bound(momentum,base[j]);
    }
    Real result=0;
    for(std::size_t i=0;i<base.size();++i) for(std::size_t j=i;j<base.size();++j) {
      const Real u=(i+j)*h;
      Real fraction=(u>mass_lower && u<mass_upper) ? 1 : 0;
      if(average_mass_cells && i+j>0)
        fraction=std::clamp((mass_upper-u+h/2)/h,Real{0},Real{1})
            -std::clamp((mass_lower-u+h/2)/h,Real{0},Real{1});
      const Real angular=(std::min(tight[i],base[j])+std::min(base[i],tight[j]))/2;
      const Real weight=first_weights[i]*second_weights[j]
          +(i==j ? 0 : first_weights[j]*second_weights[i]);
      result+=weight*angular*fraction;
    }
    return result;
  }
  std::vector<Real> first_prefix{0}, second_prefix{0}, losses;
  const auto cumulative = [&](Real cut,const std::vector<Real>& weights,
                              const std::vector<Real>& prefix) {
    if(std::isinf(cut)) return cut>0 ? prefix.back() : Real{0};
    const Real lo=cut-h/2,hi=cut+h/2;
    const auto begin=static_cast<std::size_t>(std::upper_bound(losses.begin(),losses.end(),lo)-losses.begin());
    const auto end=static_cast<std::size_t>(std::lower_bound(losses.begin(),losses.end(),hi)-losses.begin());
    Real value=prefix[begin];
    // Only one or two nodes can lie in this transition on a uniform grid.
    for(std::size_t k=begin;k<end;++k)
      value+=weights[k]*std::clamp((hi-losses[k])/h,Real{0},Real{1});
    return value;
  };
  Real result = 0;
  for (std::size_t j = 0; j < first_weights.size(); ++j) {
    const Real u = j*h;
    const Real log_ratio = cuts.pt_min == 0 ? -std::numeric_limits<Real>::infinity()
        : std::log(2*cuts.pt_min/mass)+u;
    const Real ratio = std::exp(log_ratio);
    if (ratio >= 1) break; // smaller z cannot pass both pt cuts
    const Real angular = std::min(eta_bound, std::sqrt((1-ratio)*(1+ratio)));
    // i<j => z_i>z_j; the angular cut depends only on z_j. The mass window
    // selects an interval in u_i+u_j. Sum both leg assignments bilinearly.
    if(average_mass_cells) {
      result+=angular*(first_weights[j]*(cumulative(mass_upper-u,second_weights,second_prefix)
          -cumulative(mass_lower-u,second_weights,second_prefix))
          +second_weights[j]*(cumulative(mass_upper-u,first_weights,first_prefix)
          -cumulative(mass_lower-u,first_weights,first_prefix)));
      if(j>0) {
        const Real cell=std::clamp((mass_upper-2*u+h/2)/h,Real{0},Real{1})
            -std::clamp((mass_lower-2*u+h/2)/h,Real{0},Real{1});
        result+=first_weights[j]*second_weights[j]*angular*cell;
      } else if(0>mass_lower && 0<mass_upper)
        result+=first_weights[j]*second_weights[j]*angular;
    } else {
      const auto first = std::upper_bound(losses.begin(), losses.end(), mass_lower-u);
      const auto last = std::lower_bound(losses.begin(), losses.end(), mass_upper-u);
      const auto begin_index = static_cast<std::size_t>(first-losses.begin());
      const auto end_index = static_cast<std::size_t>(last-losses.begin());
      result += angular*(
          first_weights[j]*(second_prefix[end_index]-second_prefix[begin_index])+
          second_weights[j]*(first_prefix[end_index]-first_prefix[begin_index]));
      if (2*u > mass_lower && 2*u < mass_upper)
        result += first_weights[j]*second_weights[j]*angular;
    }
    losses.push_back(u);
    first_prefix.push_back(first_prefix.back()+first_weights[j]);
    second_prefix.push_back(second_prefix.back()+second_weights[j]);
  }
  return result;
}

FFDecayMatch decay_ff_ll_match(Real mass, const FiducialCuts& cuts,
                              const LLBoundaryNonSinglet& evolution,
                              const IntegrationOptions& options) {
  const auto& ll = evolution.evolution();
  const Real a = evolution.parameters().expansion_parameter();
  FFDecayMatch result{};
  result.born = isotropic_acceptance(mass,1,1,cuts);
  result.ll = decay_ll_acceptance(mass,cuts,ll);
  result.evolved_boundary_coefficient = 2*decay_lattice_acceptance(
      mass,cuts,ll.grid().spacing(),ll.weights(),evolution.evolved_boundary_weights());
  result.resummed_with_boundary = result.ll+a*result.evolved_boundary_coefficient;
  result.lattice_overlap_coefficient = 2*evolution.first_order_apply(
      [&](Real z) { return isotropic_acceptance(mass,z,1,cuts); });
  result.fixed_order_coefficient = decay_ff_coefficient_order1(
      mass,cuts,evolution.parameters(),options);
  result.additive_correction = a*(result.fixed_order_coefficient-result.lattice_overlap_coefficient);
  result.matched = result.resummed_with_boundary+result.additive_correction;
  return result;
}
} // namespace leptonff
