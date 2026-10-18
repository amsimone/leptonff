#include "leptonff/hard_dressing.hpp"
#include "leptonff/vector_decay.hpp"
#include "cut_breakpoints.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace leptonff {
namespace {
struct Sum {
  Real value = 0, correction = 0;
  void add(Real x) {
    const Real y = x-correction, next = value+y;
    correction = (next-value)-y;
    value = next;
  }
};

// Primitive of the hard mass-spectrum kernel times dz, evaluated at z=e^-u.
// Stable logarithms cover both ends of the loss lattice.
Real hard_primitive(Real u, Real log_scale) {
  const Real z = std::exp(-u), t = -std::expm1(-u);
  const Real logt = u < .5L ? std::log(t) : std::log1p(-z);
  const Real p_primitive = -2*logt-z-z*z/2;
  if (u >= 1) {
    // Subtract the z=0 primitive analytically. Otherwise tiny far-tail cells
    // become differences of O(1) constants and contaminate tilted wrap bounds.
    return 2*(log_scale*p_primitive+(2*logt+z+z*z/2)*u
        -2*real_dilogarithm(z)+4*z+3*z*z/4-2*logt*logt
        +(3-2*z-z*z)*logt);
  }
  const Real logz_primitive = 2*real_dilogarithm(t)
      +(z+z*z/2)*u+z+z*z/4;
  const Real logsoft_primitive = -2*logt*logt
      +(4*t-t*t)*logt-4*t+t*t/2;
  const Real zeta2 = std::pow(std::acos(-1.0L),2)/6;
  return 2*(log_scale*p_primitive+logz_primitive+logsoft_primitive)-4*zeta2+7;
}

void compatible(Real mass, const FiducialCuts& cuts, const QEDParameters& q,
                 const FiducialVectorOrder1& fixed) {
  const auto& p = fixed.parameters;
  if (mass != fixed.parent_mass || cuts.pt_min != fixed.cuts.pt_min ||
      cuts.eta_max != fixed.cuts.eta_max || cuts.mass_min != fixed.cuts.mass_min ||
      cuts.mass_max != fixed.cuts.mass_max || q.alpha != p.alpha ||
      q.muon_mass != p.muon_mass || q.mu_0 != p.mu_0 || q.mu_frag != p.mu_frag)
    throw std::invalid_argument("Hard coefficient and LL dressing configurations differ");
}
Real pair_geometry(Real mass, const FiducialCuts& cuts, Real f1, Real f2,
                    const HardDecayOptions& options,const SingleLeptonSelection& selection = {}) {
  if (cuts.pt_min == 0 && std::tanh(cuts.eta_max) == 1 && selection.cut==SingleLeptonCut::none) return 0;
  const Real product = f1*f2;
  const Real pt_ratio = 2*cuts.pt_min/mass;
  if (product == 0 || f1 <= pt_ratio || f2 <= pt_ratio) return 0;
  const Real lower = std::pow(cuts.mass_min/mass,2)/product;
  const Real upper = std::min(1-options.soft_cutoff,
                              std::pow(cuts.mass_max/mass,2)/product);
  if (lower >= upper) return 0;
  std::vector<Real> outer{pt_ratio/f1,pt_ratio/f2,
                          pt_ratio*(1/f1+1/f2)-1};
  Real eta_energy = 0;
  if (cuts.pt_min > 0 && std::tanh(cuts.eta_max) < 1) {
    eta_energy = cuts.pt_min*std::cosh(cuts.eta_max);
    outer.push_back(2*eta_energy/(mass*f1));
    outer.push_back(2*eta_energy/(mass*f2));
  }
  std::vector<Real> energies{cuts.pt_min,eta_energy};
  if(selection.cut==SingleLeptonCut::pt_above) {
    energies.push_back(selection.threshold);
    if(std::tanh(cuts.eta_max)<1) energies.push_back(selection.threshold*std::cosh(cuts.eta_max));
    outer.push_back(2*(cuts.pt_min+selection.threshold)/mass-1);
  }
  if(selection.cut==SingleLeptonCut::absolute_eta_below && std::tanh(selection.threshold)<1)
    energies.push_back(cuts.pt_min*std::cosh(selection.threshold));
  if(selection.cut!=SingleLeptonCut::none)
    for(std::size_t k=2;k<energies.size();++k)
      for(Real f:{f1,f2}) {
        outer.push_back(2*energies[k]/(mass*f));
        outer.push_back(4*energies[k]/(mass*f)-1);
      }
  if(selection.cut==SingleLeptonCut::pt_above && selection.threshold>0) {
    const Real r=cuts.pt_min/selection.threshold;
    for(Real ratio:{f1/f2,f2/f1}) {outer.push_back(r*ratio);outer.push_back(ratio/r);}
  }
  return integrate([&](Real z) {
    const Real b12 = isotropic_acceptance(mass,z*f1,f2,cuts,selection);
    const Real b21 = isotropic_acceptance(mass,z*f2,f1,cuts,selection);
    std::vector<Real> inner;
    for (Real f : {f1,f2})
      for (Real energy : energies) {
        if (energy == 0) continue;
        const Real c = std::abs(((1+z)-4*energy/(mass*f))/(1-z));
        if (c > 0 && c < 1) inner.push_back(std::sqrt(1-c));
      }
    return integrate([&](Real u) {
      const Real c = 1-u*u;
      const Real a12 = vector_fragmented_acceptance(mass,z,c,f1,f2,cuts,selection);
      const Real a21 = f1 == f2 ? a12
          : vector_fragmented_acceptance(mass,z,c,f2,f1,cuts,selection);
      return 2*u*vector_scalar_real_density(mass,0,z,c)*(a12-b12+a21-b21);
    },std::sqrt(options.collinear_cutoff),1,options.integration,
        detail::resolved_cut_breakpoints(inner,std::sqrt(options.collinear_cutoff),1)).value;
  },lower,upper,options.integration,detail::resolved_cut_breakpoints(outer,lower,upper)).value;
}

} // namespace

Real vector_muon_pair_geometry(Real mass, const FiducialCuts& cuts,
    Real f1, Real f2, const HardDecayOptions& options,const SingleLeptonSelection& selection) {
  cuts.validate(); options.validate();
  selection.validate();
  if (!std::isfinite(mass) || mass <= 0 || !std::isfinite(f1) ||
      !std::isfinite(f2) || f1 < 0 || f1 > 1 || f2 < 0 || f2 > 1)
    throw std::invalid_argument("Invalid finite muon-pair geometry point");
  return pair_geometry(mass,cuts,f1,f2,options,selection);
}

LLQuadrature compress_ll_quadrature(const LLNonSinglet& ll, std::size_t panels) {
  if (panels < 1 || panels > 256)
    throw std::invalid_argument("LL quadrature panels must be in [1,256]");
  LLQuadrature result;
  result.nodes.push_back({1,ll.weights()[0]});
  Sum positive, negative, negative_z;
  for (std::size_t j = 1; j < ll.weights().size(); ++j) {
    const Real w = ll.weights()[j];
    if (w > 0) positive.add(w);
    else if (w < 0) { negative.add(w); negative_z.add(w*ll.z(j)); }
  }
  result.negative_weight_mass = -negative.value+std::max(Real{0},-ll.weights()[0]);
  if (positive.value > 0) {
    // Square-root loss panels resolve the endpoint while keeping nodes in the
    // finite-loss region where pT and mass thresholds make quantile rules coarse.
    std::vector<Sum> weights(panels), m1(panels), m2(panels), m3(panels);
    std::vector<Real> centres(panels);
    for (std::size_t k = 0; k < panels; ++k)
      centres[k] = 1-(std::pow(Real(k)/panels,2)+std::pow(Real(k+1)/panels,2))/2;
    for (std::size_t j = 1; j < ll.weights().size(); ++j) {
      const Real w = ll.weights()[j];
      if (w <= 0) continue;
      const Real z = ll.z(j);
      const auto k = std::min(panels-1,static_cast<std::size_t>(
          panels*std::sqrt(std::max(Real{0},1-z))));
      const Real x = z-centres[k];
      weights[k].add(w);
      m1[k].add(w*x); m2[k].add(w*x*x); m3[k].add(w*x*x*x);
    }
    for (std::size_t k = 0; k < panels; ++k) {
      const Real w = weights[k].value;
      if (w <= 0) continue;
      const Real mean = m1[k].value/w, second = m2[k].value/w;
      const Real variance = std::max(Real{0},second-mean*mean);
      if (variance < 1e-26L) {
        result.nodes.push_back({centres[k]+mean,w});
        continue;
      }
      // Two-node Gauss rule for the positive discrete measure in this panel.
      // Centred moments keep its Jacobi coefficients stable near z=1.
      const Real a1 = (m3[k].value/w-2*mean*second+mean*mean*mean)/variance;
      const Real gap = std::hypot(mean-a1,2*std::sqrt(variance));
      const Real lower = (mean+a1-gap)/2, upper = (mean+a1+gap)/2;
      const Real weight_lower = w*(upper-mean)/gap;
      result.nodes.push_back({centres[k]+upper,w-weight_lower});
      result.nodes.push_back({centres[k]+lower,weight_lower});
    }
  }
  if (negative.value < 0)
    result.nodes.push_back({negative_z.value/negative.value,negative.value});
  Sum number, moment, second_moment, third_moment;
  for (const auto& node : result.nodes) {
    if (!std::isfinite(node.fraction) || node.fraction < 0 || node.fraction > 1 ||
        !std::isfinite(node.weight))
      throw std::runtime_error("Invalid compressed LL node");
    number.add(node.weight);
    moment.add(node.weight*node.fraction);
    second_moment.add(node.weight*node.fraction*node.fraction);
    third_moment.add(node.weight*std::pow(node.fraction,3));
  }
  result.number_residual = number.value-ll.moment(1);
  result.first_moment_residual = moment.value-ll.moment(2);
  result.second_moment_residual = second_moment.value-ll.moment(3);
  result.third_moment_residual = third_moment.value-ll.moment(4);
  if (std::abs(result.number_residual) > 1e-12L ||
      std::abs(result.first_moment_residual) > 1e-12L ||
      std::abs(result.second_moment_residual) > 1e-12L ||
      std::abs(result.third_moment_residual) > 1e-12L)
    throw std::runtime_error("LL compression failed number/first-moment conservation");
  return result;
}

void HardDressingOptions::validate() const {
  if (ff_panels < 1 || ff_panels > 256)
    throw std::invalid_argument("Hard dressing FF panels must be in [1,256]");
  geometry.validate();
}

HardProjection evolve_hard_projection(Real mass, Real mu_frag,
                                      const LLNonSinglet& ll) {
  if (!std::isfinite(mass) || mass <= 0 || !std::isfinite(mu_frag) || mu_frag <= 0)
    throw std::invalid_argument("Hard projection requires positive finite scales");
  HardProjection result;
  const auto& grid = ll.grid();
  const Real h = grid.spacing(), log_scale = 2*(std::log(mass)-std::log(mu_frag));
  result.initial_weights.resize(grid.bins,0);
  Sum rate, tilted;
  for (std::size_t j = 1; j < grid.bins; ++j) {
    const Real lo = std::max(grid.soft_log_cutoff,(j-.5L)*h);
    const Real hi = (j+.5L)*h;
    if (lo >= hi) continue;
    const Real w = hard_primitive(lo,log_scale)-hard_primitive(hi,log_scale);
    if (!std::isfinite(w)) throw std::runtime_error("Non-finite hard lattice cell");
    result.initial_weights[j] = w;
    rate.add(w);
    tilted.add(std::abs(w)*std::exp(.9L*j*h));
  }
  result.initial_weights[0] = 1.5L-rate.value;
  tilted.add(std::abs(result.initial_weights[0]));
  const LLNonSinglet pair_ll(2*ll.time(),grid);
  result.wrap_coefficient_bound = pair_ll.diagnostics().wrap_probability_bound*tilted.value;
  const Real end = (grid.bins-.5L)*h, e = std::exp(-end);
  result.omitted_tail_coefficient_bound = 2*(1+e*e)/(1-e)*e
      *(std::abs(log_scale)+end+1+2*e/(1-e));
  auto convolved = circular_lattice_convolution(result.initial_weights,ll.weights());
  result.evolved_weights = std::move(convolved.weights);
  result.imaginary_weight_mass = convolved.imaginary_weight_mass;
  Sum number;
  for (Real w : result.evolved_weights) number.add(w);
  result.number_residual = number.value-1.5L;
  if (std::abs(result.number_residual) > 1e-7L)
    throw std::runtime_error("Hard projection number residual exceeds tolerance");
  return result;
}

HardDressedVectorMatch fiducial_vector_hard_dressed_match(Real mass,
    const FiducialCuts& cuts, const LLBoundaryNonSinglet& evolution,
    const FiducialVectorOrder1& fixed, const HardDressingOptions& options,
    const IntegrationOptions& distribution_options) {
  compatible(mass,cuts,evolution.parameters(),fixed);
  options.validate();
  const auto& ll = evolution.evolution();
  const Real a = evolution.parameters().expansion_parameter();
  HardDressedVectorMatch result{};
  result.projection = evolve_hard_projection(mass,evolution.parameters().mu_frag,ll);
  if (a*(result.projection.wrap_coefficient_bound
      +2*result.projection.omitted_tail_coefficient_bound) > ll.grid().max_wrap_probability)
    throw std::invalid_argument("Hard correction wrap/tail bound too large: increase log_span");
  result.quadrature = compress_ll_quadrature(ll,options.ff_panels);
  const auto& p = result.projection;
  std::vector<Real> identity(ll.grid().bins,0);
  identity[0] = 1;
  const Real bare_projection = decay_lattice_acceptance(mass,cuts,ll.grid().spacing(),
                                                        p.initial_weights,identity);
  const Real dressed_projection = decay_lattice_acceptance(mass,cuts,ll.grid().spacing(),
                                                           p.evolved_weights,ll.weights());
  result.projected_dressing_shift = dressed_projection-bare_projection;
  const auto& nodes = result.quadrature.nodes;
  Sum geometry_shift;
  for (std::size_t i = 0; i < nodes.size(); ++i)
    for (std::size_t j = i; j < nodes.size(); ++j) {
      Real weight = nodes[i].weight*nodes[j].weight*(i == j ? 1 : 2);
      if (i == 0 && j == 0) weight -= 1; // correlated zero-evolution subtraction
      if (weight == 0) continue;
      geometry_shift.add(weight*pair_geometry(mass,cuts,nodes[i].fraction,
                                               nodes[j].fraction,options.geometry));
    }
  result.geometry_dressing_shift = geometry_shift.value;
  result.undressed_hard_coefficient = fixed.hard_coefficient;
  result.dressed_hard_coefficient = fixed.hard_coefficient
      +result.projected_dressing_shift+result.geometry_dressing_shift;
  const auto ff = decay_ff_ll_match(mass,cuts,evolution,distribution_options);
  result.ff_matched = ff.matched;
  result.additive_matched = ff.matched+a*fixed.hard_coefficient;
  result.matched = ff.matched+a*result.dressed_hard_coefficient;
  return result;
}
} // namespace leptonff
