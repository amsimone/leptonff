#include "flavour_evolution.hpp"
#include "fft.hpp"
#include <complex>

namespace leptonff { namespace detail {
namespace {
struct ProjectedBasis {
  std::vector<Real> weights;
  std::vector<Complex> frequency;
  Real tail=0,resolved_norm=0,weighted_resolved_norm=0;
};
// Same energy-weighted linear-in-z deposition as the legacy backend. Keep
// this independent implementation auditable against moment and legacy tests.
ProjectedBasis project_basis(const Distribution& d,const QEDXSpaceGrid& grid,const IntegrationOptions& options) {
  const auto n=grid.bins;const Real h=grid.spacing(),z1=std::exp(-h),cap=std::exp(-h*(n-1));
  ProjectedBasis p;p.weights.resize(n);
  p.weights[1]=d.apply([=](Real z){return z*(1-z)/(1-z1);},z1,options);
  for(std::size_t j=1;j+1<n;++j) {
    const Real za=std::exp(-h*j),zb=std::exp(-h*(j+1));
    const auto density=[&](Real u){const Real z=std::exp(-u);return z*z*d.resolved_density(z);};
    const Real total=integrate(density,j*h,(j+1)*h,options).value;
    const Real right=integrate([&](Real u){return density(u)*(za-std::exp(-u))/(za-zb);},j*h,(j+1)*h,options).value;
    p.weights[j]+=total-right;p.weights[j+1]+=right;
  }
  auto tail_options=options;tail_options.absolute_tolerance*=cap;
  const Real tail=integrate([&](Real t){const Real z=cap*t;return cap*z*d.resolved_density(z);},0,1,tail_options).value;
  const auto absolute=integrate([&](Real t){const Real z=cap*t;return cap*z*std::abs(d.resolved_density(z));},0,1,tail_options);
  p.tail=absolute.value+absolute.estimated_error;
  Real resolved=0;for(std::size_t k=1;k<n;++k) {
    resolved+=p.weights[k];p.resolved_norm+=std::abs(p.weights[k]);
    p.weighted_resolved_norm+=std::abs(p.weights[k])*std::exp(grid.wrap_exponent*h*k);
  }
  p.weights[0]=d.moment(2,options)-resolved-tail;
  p.frequency.assign(p.weights.begin(),p.weights.end());fft(p.frequency,false);return p;
}
using Basis=std::array<ProjectedBasis,flavour_basis_count>;
template<class T> T evaluate(const FlavourExpansion& c,const std::array<T,flavour_basis_count>& b) {
  T value=0;for(std::size_t k=0;k<c.size();++k)if(c[k]!=0)value+=c[k]*b[k];return value;
}
Real atom(const FlavourExpansion& c,const Basis& b) {
  Real value=0;for(std::size_t k=0;k<c.size();++k)value+=c[k]*b[k].weights[0];return value;
}
Real resolved_norm(const FlavourExpansion& c,const Basis& b,bool weighted) {
  Real value=0;for(std::size_t k=0;k<c.size();++k)value+=std::abs(c[k])*(weighted?b[k].weighted_resolved_norm:b[k].resolved_norm);return value;
}
Real tail_norm(const FlavourExpansion& c,const Basis& b) {
  Real value=0;for(std::size_t k=0;k<c.size();++k)value+=std::abs(c[k])*b[k].tail;return value;
}
Real log_norm(const std::vector<FlavourExpansion>& matrix,std::size_t n,const Basis& b,bool weighted) {
  Real maximum=-1e300L;
  for(std::size_t i=0;i<n;++i) {
    Real row=0;
    for(std::size_t j=0;j<n;++j) {
      const auto& c=matrix[i*n+j];row+=(i==j?atom(c,b):std::abs(atom(c,b)))+resolved_norm(c,b,weighted);
    }
    maximum=std::max(maximum,row);
  }
  return maximum;
}
Real matrix_tail(const std::vector<FlavourExpansion>& matrix,std::size_t n,const Basis& b) {
  Real maximum=0;for(std::size_t i=0;i<n;++i) {
    Real row=0;for(std::size_t j=0;j<n;++j)row+=tail_norm(matrix[i*n+j],b);maximum=std::max(maximum,row);
  }
  return maximum;
}
struct PhaseExpansion {
  Real length=0,beta0=0,beta1=0,photon_matching=0,strong_initial=0,strong_beta=0;
  std::vector<FlavourExpansion> p0,p1,qcd0;
  std::vector<std::pair<std::size_t,FlavourExpansion>> heavy;
};
using SelectedSpectrum=std::array<std::array<std::vector<Complex>,3>,3>;
void put(std::vector<Complex>& v,std::size_t k,Complex z,std::size_t n) {
  v[k]=z;if(k && k!=n/2)v[n-k]=std::conj(z);
}
void inverse_selected(SelectedSpectrum& f,LatticeMatrix& result) {
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
    fft(f[i][j],true);result[i][j].resize(f[i][j].size());
    for(std::size_t k=0;k<f[i][j].size();++k)result[i][j][k]=f[i][j][k].real();
  }
}
} // namespace
FlavourXSpaceResult evolve_flavour_xspace(const MomentEvolutionParameters& p,const QEDXSpaceGrid& grid,
    const MomentEvolutionOptions& options) {
  const auto& s=*p.flavour_scheme;const auto n=s.size(),bins=grid.bins;const auto selected=s.selected();
  const Real a=p.alpha_msbar_at_mu0/(2*std::acos(-1.L));
  Basis basis;for(std::size_t k=0;k<basis.size();++k)basis[k]=project_basis(flavour_basis_distribution(k),grid,options.kernel_integration);
  std::vector<std::array<FlavourExpansion,3>> boundary(n);
  const auto initial=s.active(p.mu_0),final=s.active(p.mu_frag);
  for(std::size_t i=0;i<n;++i)for(unsigned j=0;j<3;++j)
    boundary[i][j]=flavour_boundary_expansion(s,initial,i,selected[j],p.mu_0);
  Real boundary_norm0=0,boundary_norm_sigma=0,boundary_tail=0;
  for(std::size_t i=0;i<n;++i) {
    Real norm0=0,norm_sigma=0,tail=0;
    for(unsigned j=0;j<3;++j) {
      const Real endpoint=std::abs((i==selected[j]?1:0)+a*atom(boundary[i][j],basis));
      norm0+=endpoint+a*resolved_norm(boundary[i][j],basis,false);
      norm_sigma+=endpoint+a*resolved_norm(boundary[i][j],basis,true);
      tail+=tail_norm(boundary[i][j],basis);
    }
    boundary_norm0=std::max(boundary_norm0,norm0);boundary_norm_sigma=std::max(boundary_norm_sigma,norm_sigma);boundary_tail=std::max(boundary_tail,tail);
  }
  std::vector<PhaseExpansion> expansions;
  Real log_growth=0,log_growth_sigma=0,missing=0;
  for(const auto& interval:flavour_intervals(p)) {
    PhaseExpansion phase;phase.length=2*std::log(interval.upper/interval.lower);
    const auto sums=qed_charge_sums(s,interval.active);phase.beta0=sums.beta0();phase.beta1=sums.beta1();
    phase.p0.resize(n*n);phase.p1.resize(n*n);
    if(p.qcd){phase.qcd0.resize(n*n);phase.strong_initial=run_qcd_coupling(p,interval.lower);
      phase.strong_beta=11.L/2;for(std::size_t f=0;f<interval.active.size();++f)if(interval.active[f] && s.fermions[f].colours==3)phase.strong_beta-=1.L/3;}
    for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<n;++j) {
      phase.p0[i*n+j]=flavour_kernel_expansion(s,interval.active,i,j,0);
      if(p.qcd)phase.qcd0[i*n+j]=qcd_kernel_expansion(s,interval.active,i,j);
      phase.p1[i*n+j]=flavour_kernel_expansion(s,interval.active,i,j,1);
    }
    // Running integrals in the LOW theory, without a threshold at this edge.
    auto low=p;low.mu_0=interval.lower;low.mu_frag=interval.upper;
    low.alpha_msbar_at_mu0=2*std::acos(-1.L)*run_flavour_coupling(p,interval.lower).a_final;
    // An endpoint threshold changes a but not these open-interval integrals.
    const auto flow=run_flavour_coupling(low,interval.upper);
    log_growth+=flow.integral_a*log_norm(phase.p0,n,basis,false)+flow.integral_a2*log_norm(phase.p1,n,basis,false);
    log_growth_sigma+=flow.integral_a*log_norm(phase.p0,n,basis,true)+flow.integral_a2*log_norm(phase.p1,n,basis,true);
    missing+=flow.integral_a*matrix_tail(phase.p0,n,basis)+flow.integral_a2*matrix_tail(phase.p1,n,basis);
    if(p.qcd) {
      const Real I=integral_qcd_coupling(p,interval.lower,interval.upper);
      log_growth+=I*log_norm(phase.qcd0,n,basis,false);
      log_growth_sigma+=I*log_norm(phase.qcd0,n,basis,true);
      missing+=I*matrix_tail(phase.qcd0,n,basis);
    }
    Real match_norm0=1,match_norm_sigma=1,match_tail=0;
    // flow.a_final is ABOVE matching; reconstruct the below-threshold a.
    for(auto f:interval.activating)phase.photon_matching-=flavour_matching_expansion(s,f,s.photon(),s.photon(),interval.upper)[Delta];
    if(p.hadronic_matching && interval.upper==p.hadronic_matching->scale)phase.photon_matching+=p.hadronic_matching->coefficient;
    const Real below=p.coupling==QEDCouplingMode::fixed?flow.a_final:flow.a_final/(1+flow.a_final*phase.photon_matching);
    match_norm0=std::max(match_norm0,std::abs(1-below*phase.photon_matching));match_norm_sigma=match_norm0;
    for(auto f:interval.activating) {
      auto c=flavour_matching_expansion(s,f,2*f,s.photon(),interval.upper);phase.heavy.push_back({2*f,c});
      match_norm0=std::max(match_norm0,below*(std::abs(atom(c,basis))+resolved_norm(c,basis,false)));
      match_norm_sigma=std::max(match_norm_sigma,below*(std::abs(atom(c,basis))+resolved_norm(c,basis,true)));
      match_tail=std::max(match_tail,below*tail_norm(c,basis));
    }
    log_growth+=std::log(match_norm0);log_growth_sigma+=std::log(match_norm_sigma);missing+=match_tail;
    expansions.push_back(std::move(phase));
  }
  FlavourXSpaceResult r;
  r.diagnostics.weighted_wrap_bound=std::exp(-grid.wrap_exponent*grid.log_span+log_growth_sigma)*boundary_norm_sigma;
  r.diagnostics.weighted_tail_bound=std::exp(log_growth)*(boundary_norm0*std::expm1(missing)+a*boundary_tail*std::exp(missing));
  if(!std::isfinite(r.diagnostics.weighted_wrap_bound) || !std::isfinite(r.diagnostics.weighted_tail_bound) ||
      r.diagnostics.weighted_wrap_bound>grid.max_wrap_bound || r.diagnostics.weighted_tail_bound>grid.max_tail_bound)
    throw std::invalid_argument("Variable-flavour lattice fails weighted wrap/tail bounds");
  std::vector<std::array<std::vector<Complex>,3>> fd(n);
  for(auto& row:fd)for(auto& v:row)v.resize(bins);
  SelectedSpectrum fe{},fc1{},fc2{};
  for(auto* f:{&fe,&fc1,&fc2})for(auto& row:*f)for(auto& v:row)v.resize(bins);
  const auto combine_weights=[&](const FlavourExpansion& c) {
    std::vector<Real> w(bins);for(std::size_t k=0;k<c.size();++k)if(c[k]!=0)
      for(std::size_t j=0;j<bins;++j)w[j]+=c[k]*basis[k].weights[j];
    return w;
  };
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
    r.boundary1[i][j]=combine_weights(boundary[selected[i]][j]);
    r.kernel0[i][j]=combine_weights(flavour_kernel_expansion(s,final,selected[i],selected[j],0));
    r.kernel1[i][j]=combine_weights(flavour_kernel_expansion(s,final,selected[i],selected[j],1));
  }
  for(std::size_t k=0;k<=bins/2;++k) {
    std::array<Complex,flavour_basis_count> frequency;
    for(std::size_t b=0;b<frequency.size();++b)frequency[b]=basis[b].frequency[k];
    FlavourColumns<Complex> b1(n);
    for(std::size_t i=0;i<n;++i)for(unsigned j=0;j<3;++j)b1[i][j]=evaluate(boundary[i][j],frequency);
    std::vector<FlavourPhase<Complex>> phases;
    for(const auto& expansion:expansions) {
      FlavourPhase<Complex> phase;phase.length=expansion.length;phase.beta0=expansion.beta0;phase.beta1=expansion.beta1;
      phase.photon_matching=expansion.photon_matching;phase.strong_initial=expansion.strong_initial;phase.strong_beta=expansion.strong_beta;
      if(!expansion.qcd0.empty()){phase.qcd0.resize(n*n);for(std::size_t i=0;i<n*n;++i)phase.qcd0[i]=evaluate(expansion.qcd0[i],frequency);}
      phase.p0.resize(n*n);phase.p1.resize(n*n);
      for(std::size_t i=0;i<n*n;++i){phase.p0[i]=evaluate(expansion.p0[i],frequency);phase.p1[i]=evaluate(expansion.p1[i],frequency);}
      for(const auto& h:expansion.heavy)phase.heavy_matching.push_back({h.first,evaluate(h.second,frequency)});
      phases.push_back(std::move(phase));
    }
    std::size_t steps;Real shift;
    const auto state=evolve_flavour_columns(s,a,b1,phases,p.coupling,options,steps,shift);
    r.diagnostics.maximum_steps=std::max(r.diagnostics.maximum_steps,steps);
    r.diagnostics.maximum_frequency_step_shift=std::max(r.diagnostics.maximum_frequency_step_shift,shift);
    for(std::size_t i=0;i<n;++i)for(unsigned j=0;j<3;++j)put(fd[i][j],k,state.d[i][j],bins);
    for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
      put(fe[i][j],k,state.e[selected[i]][j],bins);put(fc1[i][j],k,state.c1[selected[i]][j],bins);put(fc2[i][j],k,state.c2[selected[i]][j],bins);
    }
  }
  r.full_weights.resize(n);
  for(std::size_t i=0;i<n;++i)for(unsigned j=0;j<3;++j) {
    fft(fd[i][j],true);auto& w=r.full_weights[i][j];w.resize(bins);
    for(std::size_t k=0;k<bins;++k) {
      w[k]=fd[i][j][k].real();r.diagnostics.imaginary_mass+=std::abs(fd[i][j][k].imag());
      if(!std::isfinite(w[k]))throw std::runtime_error("Nonfinite variable-flavour lattice weight");
      if(w[k]<0)r.diagnostics.negative_weight_mass-=w[k];
    }
  }
  if(r.diagnostics.imaginary_mass>1e-8L)throw std::runtime_error("Variable-flavour inverse FFT imaginary residual too large");
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)r.weights[i][j]=r.full_weights[selected[i]][j];
  inverse_selected(fe,r.propagator);inverse_selected(fc1,r.c1);inverse_selected(fc2,r.c2);
  r.diagnostics.weighted_step_shift_bound=3*n*std::sqrt(Real(bins))*r.diagnostics.maximum_frequency_step_shift;
  return r;
}
}} // namespace leptonff::detail
