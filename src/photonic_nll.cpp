#include "leptonff/photonic_nll.hpp"
#include "leptonff/nlo_lepton.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace leptonff {
Real photonic_source_number_moment() {
  constexpr Real zeta3=1.20205690315959428539973816151145L;
  const Real pi=std::acos(-1.L);
  return -13.L/4+pi*pi/2-2*zeta3;
}
Distribution unitary_photonic_nlo_kernel() {
  // E_gamma(N)/E_gamma(1) to NLL, not a change to the MSbar source G.
  return nlo_lepton_component(NLOLeptonComponent::photonic)
    +Distribution(-photonic_source_number_moment());
}
PhotonicNLL::PhotonicNLL(const QEDParameters& p,const LLGrid& g)
    :PhotonicNLL(LLBoundaryNonSinglet(p,g)) {}
PhotonicNLL::PhotonicNLL(const LLBoundaryNonSinglet& baseline):baseline_(baseline) {
  const auto& u=baseline_.evolution();const auto& g=u.grid();
  const auto& p=baseline_.parameters();
  const Real a=p.expansion_parameter(),L=2*std::log(p.mu_frag/p.mu_0),h=g.spacing();
  const auto source=nlo_lepton_component(NLOLeptonComponent::photonic);
  const auto density=[&](Real v) {const Real z=std::exp(-v);return z*source.resolved_density(z);};
  kernel_.resize(g.bins,0);
  Real sum=0,tilted=0;
  const IntegrationOptions options{1e-14L,2e-12L,8192};
  for(std::size_t k=1;k<g.bins;++k) {
    const Real lo=std::max(g.soft_log_cutoff,(k-.5L)*h),hi=(k+.5L)*h;
    if(lo>=hi)continue;
    kernel_[k]=integrate(density,lo,hi,options).value;
    sum+=kernel_[k];tilted+=std::abs(kernel_[k])*std::exp(.9L*k*h);
  }
  // For G-G(1)*delta, endpoint = -integral G_regular. Unresolved cells
  // are assigned to the atom, exactly as for P0; refinement controls this.
  kernel_[0]=-sum;tilted+=std::abs(kernel_[0]);
  const auto convolution=circular_lattice_convolution(u.weights(),kernel_);
  evolved_=convolution.weights;
  diagnostics_.source_number_moment=photonic_source_number_moment();
  diagnostics_.raw_inclusive_factor=1+a*a*L*diagnostics_.source_number_moment;
  diagnostics_.kernel_number_residual=std::accumulate(kernel_.begin(),kernel_.end(),Real{0});
  diagnostics_.evolved_number_residual=std::accumulate(evolved_.begin(),evolved_.end(),Real{0});
  diagnostics_.imaginary_weight_mass=convolution.imaginary_weight_mass;
  diagnostics_.wrap_bound=a*a*L*u.diagnostics().wrap_probability_bound*tilted;
  // For v>=20 the source density is bounded conservatively by
  // 20 exp(-v) (1+v^2); two accounts for replacing its omitted tail by delta.
  const Real end=(g.bins-.5L)*h;
  if(end<20)throw std::invalid_argument("Photonic NLL tail bound requires log span >=20 plus half a cell");
  diagnostics_.omitted_tail_bound=40*a*a*L*std::exp(-end)*(end*end+2*end+3);
  if(!std::isfinite(diagnostics_.wrap_bound) || diagnostics_.wrap_bound>g.max_wrap_probability ||
     !std::isfinite(diagnostics_.omitted_tail_bound) || diagnostics_.omitted_tail_bound>g.max_wrap_probability)
    throw std::invalid_argument("Photonic NLL domain exceeds wrap/tail tolerance: increase log span");
  insertion_=baseline_.evolved_boundary_weights();
  const auto pp=circular_lattice_convolution(u.generator_weights(),u.generator_weights());
  const auto pd=circular_lattice_convolution(u.generator_weights(),baseline_.initial_boundary_weights());
  order2_.resize(g.bins);
  for(std::size_t k=0;k<g.bins;++k) {
    insertion_[k]+=a*L*evolved_[k];
    order2_[k]=L*L*pp.weights[k]/2+L*(pd.weights[k]+kernel_[k]);
  }
}
Real PhotonicNLL::apply(const Function& f) const {
  const auto& u=baseline_.evolution();
  const Real a=baseline_.parameters().expansion_parameter();
  Real value=u.apply(f);
  for(std::size_t k=0;k<insertion_.size();++k) {
    const Real sample=f(u.z(k));
    if(!std::isfinite(sample))throw std::invalid_argument("Non-finite photonic NLL measurement");
    value+=a*insertion_[k]*sample;
  }
  if(!std::isfinite(value))throw std::runtime_error("Non-finite photonic NLL action");
  return value;
}
Real PhotonicNLL::moment(Real N) const {
  if(!std::isfinite(N) || N<1)throw std::invalid_argument("Photonic NLL lattice moments require N>=1");
  return apply([N](Real z){return std::pow(z,N-1);});
}
Real PhotonicNLL::exact_moment(Real N,const IntegrationOptions& options) const {
  if(!std::isfinite(N) || N<1)throw std::invalid_argument("Photonic NLL moments require N>=1");
  const auto& p=baseline_.parameters();const Real a=p.expansion_parameter(),L=2*std::log(p.mu_frag/p.mu_0);
  const Real d=initial_condition_coefficient(Species::muon_minus,Species::muon_minus,p).moment(N,options);
  return baseline_.evolution().exact_moment(N)*(1+a*d+a*a*L*unitary_photonic_nlo_kernel().moment(N,options));
}
} // namespace leptonff
