#include "flavour_evolution.hpp"
#include <set>

namespace leptonff {
namespace detail {
std::vector<FlavourInterval> flavour_intervals(const MomentEvolutionParameters& p) {
  if(!p.flavour_scheme)throw std::invalid_argument("Flavour evolution requires a flavour scheme");
  if(p.mu_frag<p.mu_0)throw std::invalid_argument("Variable-flavour FF evolution requires mu_frag>=mu0");
  const auto& s=*p.flavour_scheme;std::set<Real> edges;
  for(const auto& f:s.fermions)if(f.threshold()>p.mu_0 && f.threshold()<=p.mu_frag)edges.insert(f.threshold());
  if(p.qcd && p.qcd->start_scale>p.mu_0 && p.qcd->start_scale<=p.mu_frag)edges.insert(p.qcd->start_scale);
  if(p.hadronic_matching && p.hadronic_matching->scale>p.mu_0 && p.hadronic_matching->scale<=p.mu_frag)edges.insert(p.hadronic_matching->scale);
  edges.insert(p.mu_frag);Real previous=p.mu_0;std::vector<FlavourInterval> out;
  for(Real edge:edges) {
    if(edge==previous)continue;
    FlavourInterval interval{previous,edge,s.active(std::sqrt(previous*edge)),{}};
    for(std::size_t f=0;f<s.fermions.size();++f)if(s.fermions[f].threshold()==edge)interval.activating.push_back(f);
    out.push_back(interval);previous=edge;
  }
  return out;
}
namespace {
Real evaluate(const FlavourExpansion& c,const std::array<Real,flavour_basis_count>& basis) {
  Real v=0;for(std::size_t k=0;k<c.size();++k)v+=c[k]*basis[k];return v;
}
}
QEDMomentEvolution evolve_flavour_moment(Real N,const MomentEvolutionParameters& p,const MomentEvolutionOptions& options) {
  const auto& s=*p.flavour_scheme;const auto n=s.size();const auto selected=s.selected();
  std::array<Real,flavour_basis_count> basis{};
  for(std::size_t k=0;k<basis.size();++k)basis[k]=flavour_basis_distribution(k).moment(N,options.kernel_integration);
  FlavourColumns<Real> boundary(n);
  const auto initial=s.active(p.mu_0),final=s.active(p.mu_frag);
  QEDMomentEvolution r;r.N=N;r.log_scale_ratio=2*std::log(p.mu_frag/p.mu_0);
  r.a_initial=p.alpha_msbar_at_mu0/(2*std::acos(-1.L));
  for(std::size_t i=0;i<n;++i)for(unsigned j=0;j<3;++j)
    boundary[i][j]=evaluate(flavour_boundary_expansion(s,initial,i,selected[j],p.mu_0),basis);
  std::vector<FlavourPhase<Real>> phases;
  for(const auto& interval:flavour_intervals(p)) {
    const auto sums=qed_charge_sums(s,interval.active);
    FlavourPhase<Real> phase;phase.length=2*std::log(interval.upper/interval.lower);
    phase.beta0=sums.beta0();phase.beta1=sums.beta1();phase.p0.resize(n*n);phase.p1.resize(n*n);
    if(p.qcd){phase.qcd0.resize(n*n);phase.strong_initial=run_qcd_coupling(p,interval.lower);
      phase.strong_beta=11.L/2;for(std::size_t f=0;f<interval.active.size();++f)if(interval.active[f] && s.fermions[f].colours==3)phase.strong_beta-=1.L/3;}
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
      phase.p0[i*n+j]=evaluate(flavour_kernel_expansion(s,interval.active,i,j,0),basis);
      if(p.qcd)phase.qcd0[i*n+j]=evaluate(qcd_kernel_expansion(s,interval.active,i,j),basis);
      phase.p1[i*n+j]=evaluate(flavour_kernel_expansion(s,interval.active,i,j,1),basis);
    }
    for(auto f:interval.activating) {
      phase.heavy_matching.push_back({2*f,evaluate(flavour_matching_expansion(s,f,2*f,s.photon(),interval.upper),basis)});
      phase.photon_matching-=flavour_matching_expansion(s,f,s.photon(),s.photon(),interval.upper)[Delta];
    }
    if(p.hadronic_matching && interval.upper==p.hadronic_matching->scale)phase.photon_matching+=p.hadronic_matching->coefficient;
    phases.push_back(std::move(phase));
  }
  const auto state=evolve_flavour_columns(s,r.a_initial,boundary,phases,p.coupling,options,r.steps,r.step_refinement_shift);
  r.a_final=state.a;r.integral_a=state.integral_a;r.integral_a2=state.integral_a2;
  r.flavour_fragmentation=state.d;r.flavour_evolution=state.e;r.flavour_coefficient1=state.c1;r.flavour_coefficient2=state.c2;
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
    const auto row=selected[i],col=selected[j];
    r.fragmentation[i][j]=state.d[row][j];r.evolution[i][j]=state.e[row][j];
    r.boundary1[i][j]=boundary[row][j];r.coefficient1[i][j]=state.c1[row][j];r.coefficient2[i][j]=state.c2[row][j];
    r.kernel0[i][j]=evaluate(flavour_kernel_expansion(s,final,row,col,0),basis);
    r.kernel1[i][j]=evaluate(flavour_kernel_expansion(s,final,row,col,1),basis);
  }
  return r;
}
} // namespace detail
namespace {
FlavourRunningResult running_interval(Real a,Real L,QEDChargeSums sums,QEDCouplingMode mode) {
  FlavourRunningResult r;r.a_final=a;
  if(mode!=QEDCouplingMode::fixed)r.beta0_integral=sums.beta0()*L;
  if(a==0)return r;
  if(mode==QEDCouplingMode::fixed || sums.squared==0){r.integral_a=a*L;r.integral_a2=a*a*L;return r;}
  const Real b=sums.beta0(),c=mode==QEDCouplingMode::two_loop?sums.beta1():0,u0=1/a;
  Real u=u0-b*L;
  if(u<=2*std::acos(-1.L))throw std::runtime_error("Flavour running outside perturbative range");
  if(c!=0)for(int k=0;k<64;++k) {
    const Real f=(u0-u)/b-c/(b*b)*std::log1p(b*(u0-u)/(b*u+c))-L;
    const Real step=f*(b*u+c)/u;u+=step;
    if(!std::isfinite(u) || u<=2*std::acos(-1.L))throw std::runtime_error("Flavour running outside perturbative range");
    if(std::abs(step)<2e-15L*u)break;
    if(k==63)throw std::runtime_error("Flavour two-loop flow did not converge");
  }
  r.a_final=1/u;r.beta0_integral=b*L;
  const Real change=r.a_final-a;
  const Real beta_log=c==0?0:std::log1p(c*change/(b+c*a));
  r.integral_a=(std::log1p(change/a)-beta_log)/b;
  r.integral_a2=c==0?change/b:beta_log/c;
  return r;
}
}
FlavourRunningResult run_flavour_coupling(const MomentEvolutionParameters& p,Real scale) {
  p.validate();if(!p.flavour_scheme || !std::isfinite(scale) || scale<=0)throw std::invalid_argument("Invalid variable-flavour running request");
  const auto& s=*p.flavour_scheme;const bool up=scale>=p.mu_0;
  std::set<Real> points;
  for(const auto& f:s.fermions) {
    const Real t=f.threshold();if(up?(t>p.mu_0 && t<=scale):(t<=p.mu_0 && t>scale))points.insert(t);
  }
  if(p.hadronic_matching) {const Real t=p.hadronic_matching->scale;if(up?(t>p.mu_0 && t<=scale):(t<=p.mu_0 && t>scale))points.insert(t);}
  std::vector<Real> edges(points.begin(),points.end());if(!up)std::reverse(edges.begin(),edges.end());
  edges.push_back(scale);Real previous=p.mu_0;FlavourRunningResult result;
  result.a_final=p.alpha_msbar_at_mu0/(2*std::acos(-1.L));
  for(std::size_t k=0;k<edges.size();++k) {
    const Real edge=edges[k];
    if(edge!=previous) {
      const auto flow=running_interval(result.a_final,2*std::log(edge/previous),qed_charge_sums(s,s.active(std::sqrt(edge*previous))),p.coupling);
      result.a_final=flow.a_final;result.integral_a+=flow.integral_a;result.integral_a2+=flow.integral_a2;result.beta0_integral+=flow.beta0_integral;
    }
    if(k+1<edges.size()) {
      Real match=0;for(const auto& f:s.fermions)if(f.threshold()==edge && !f.pointlike)
        match+=4*f.colours*f.charge*f.charge*std::log(edge/f.mass)/3;
      if(p.hadronic_matching && edge==p.hadronic_matching->scale)match+=p.hadronic_matching->coefficient;
      if(p.coupling!=QEDCouplingMode::fixed) {
        const Real sign=up?1:-1;result.a_final/=1-sign*result.a_final*match;result.beta0_integral+=sign*match;
      }
    }
    previous=edge;
  }
  if(!std::isfinite(result.a_final) || result.a_final<0 || result.a_final>=1/(2*std::acos(-1.L)))
    throw std::runtime_error("Matched coupling outside perturbative range");
  return result;
}
} // namespace leptonff
