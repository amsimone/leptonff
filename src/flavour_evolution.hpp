#pragma once
#include "leptonff/xspace_evolution.hpp"
#include "flavour_kernels.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace leptonff { namespace detail {
template<class T> using FlavourColumns=std::vector<std::array<T,3>>;
template<class T> struct FlavourPhase {
  Real length=0,beta0=0,beta1=0,photon_matching=0,strong_initial=0,strong_beta=0;
  std::vector<T> p0,p1,qcd0; // full parent -> intermediate matrices
  std::vector<std::pair<std::size_t,T>> heavy_matching;
};
struct FlavourInterval {
  Real lower=0,upper=0;
  std::vector<bool> active;
  std::vector<std::size_t> activating;
};
std::vector<FlavourInterval> flavour_intervals(const MomentEvolutionParameters&);
template<class T> struct FlavourState {
  Real a=0,b=0,strong=0,integral_a=0,integral_a2=0;
  FlavourColumns<T> d,e,c1,c2;
};
template<class T> FlavourColumns<T> selected_identity(const QEDFlavourScheme& s) {
  FlavourColumns<T> out(s.size());const auto selected=s.selected();
  for(unsigned j=0;j<3;++j)out[selected[j]][j]=1;
  return out;
}
template<class T> FlavourState<T> flavour_add(FlavourState<T> y,const FlavourState<T>& k,Real h) {
  y.strong+=h*k.strong;
  y.a+=h*k.a;y.b+=h*k.b;y.integral_a+=h*k.integral_a;y.integral_a2+=h*k.integral_a2;
  for(std::size_t i=0;i<y.d.size();++i)for(unsigned j=0;j<3;++j) {
    y.d[i][j]+=h*k.d[i][j];y.e[i][j]+=h*k.e[i][j];y.c1[i][j]+=h*k.c1[i][j];y.c2[i][j]+=h*k.c2[i][j];
  }
  return y;
}
template<class T> FlavourState<T> flavour_derivative(const FlavourState<T>& y,
    const FlavourPhase<T>& p,const std::array<std::size_t,3>& selected,QEDCouplingMode mode) {
  if(!std::isfinite(y.a) || y.a<0 || y.a>=1/(2*std::acos(-1.L)))throw std::runtime_error("Flavour running alpha left perturbative range");
  const auto n=y.d.size();FlavourState<T> r;
  r.d.resize(n);r.e.resize(n);r.c1.resize(n);r.c2.resize(n);
  if(mode!=QEDCouplingMode::fixed){r.a=p.beta0*y.a*y.a;r.b=p.beta0;}
  if(mode==QEDCouplingMode::two_loop)r.a+=p.beta1*y.a*y.a*y.a;
  r.strong=-p.strong_beta*y.strong*y.strong;
  if(!std::isfinite(y.strong) || y.strong<0 || y.strong>=1/(2*std::acos(-1.L)))throw std::runtime_error("QCD running alpha left perturbative range");
  r.integral_a=y.a;r.integral_a2=y.a*y.a;
  for(std::size_t i=0;i<n;++i) {
    for(std::size_t j=0;j<n;++j) {
      const auto p0=p.p0[i*n+j],p1=p.p1[i*n+j],strong=p.qcd0.empty()?T{}:y.strong*p.qcd0[i*n+j],g=y.a*p0+y.a*y.a*p1+strong;
      for(unsigned c=0;c<3;++c) {
        r.d[i][c]+=g*y.d[j][c];r.e[i][c]+=g*y.e[j][c];
        if(!p.qcd0.empty()) {
          r.c1[i][c]+=strong*y.c1[j][c];r.c2[i][c]+=strong*y.c2[j][c]+p0*y.c1[j][c];
          if(j==selected[c]){r.c1[i][c]+=p0;r.c2[i][c]+=p1+y.b*p0;}
        }
      }
    }
  }
  return r;
}
template<class T> void flavour_match(FlavourState<T>& y,const FlavourPhase<T>& p,
    std::size_t gamma,QEDCouplingMode mode) {
  if(p.heavy_matching.empty() && p.photon_matching==0)return;
  const auto d=y.d[gamma],e=y.e[gamma],c1=y.c1[gamma];const Real a=y.a,b=y.b;
  for(const auto& heavy:p.heavy_matching)for(unsigned j=0;j<3;++j) {
    const auto k=heavy.second;const Real identity=j==2?1:0;
    for(auto i:{heavy.first,heavy.first+1}) {
      y.d[i][j]=a*k*d[j];y.e[i][j]=a*k*e[j];
      y.c1[i][j]=k*identity;y.c2[i][j]=k*(c1[j]+b*identity);
    }
  }
  for(unsigned j=0;j<3;++j) {
    const Real identity=j==2?1:0;
    y.d[gamma][j]*=1-a*p.photon_matching;y.e[gamma][j]*=1-a*p.photon_matching;
    y.c1[gamma][j]-=p.photon_matching*identity;
    y.c2[gamma][j]-=p.photon_matching*(c1[j]+b*identity);
  }
  if(mode!=QEDCouplingMode::fixed) {
    y.a=a/(1-a*p.photon_matching);y.b+=p.photon_matching;
    if(!std::isfinite(y.a) || y.a<0 || y.a>=1/(2*std::acos(-1.L)))throw std::runtime_error("Threshold matching left perturbative range");
  }
}
template<class T> FlavourState<T> integrate_flavour_steps(const QEDFlavourScheme& scheme,Real a,
    const FlavourColumns<T>& boundary,const std::vector<FlavourPhase<T>>& phases,
    QEDCouplingMode mode,std::size_t steps) {
  const auto identity=selected_identity<T>(scheme);FlavourState<T> y;
  y.a=a;y.d=identity;y.e=identity;y.c1=boundary;y.c2.resize(scheme.size());
  for(std::size_t i=0;i<scheme.size();++i)for(unsigned j=0;j<3;++j)y.d[i][j]+=a*boundary[i][j];
  for(const auto& p:phases) {
    const Real h=p.length/steps,before_beta=y.b;y.strong=p.strong_initial;
    const auto before_c1=y.c1,before_c2=y.c2;
    // Coefficients obey polynomial ODEs on each interval. Integrate them
    // exactly instead of letting roundoff in a large C2 drive RK refinement.
    const auto n=scheme.size();const auto selected=scheme.selected();
    if(p.qcd0.empty())for(std::size_t i=0;i<n;++i)for(unsigned c=0;c<3;++c) {
      T p0c1=0,p0square=0;
      for(std::size_t j=0;j<n;++j) {
        p0c1+=p.p0[i*n+j]*before_c1[j][c];
        p0square+=p.p0[i*n+j]*p.p0[j*n+selected[c]];
      }
      const auto source=p.p0[i*n+selected[c]];
      y.c1[i][c]=before_c1[i][c]+p.length*source;
      y.c2[i][c]=before_c2[i][c]+p.length*(p0c1+p.p1[i*n+selected[c]]+before_beta*source)
        +p.length*p.length/2*(p0square+(mode==QEDCouplingMode::fixed?0:p.beta0)*source);
    }
    for(std::size_t k=0;k<steps;++k) {
      const auto k1=flavour_derivative(y,p,scheme.selected(),mode);
      const auto k2=flavour_derivative(flavour_add(y,k1,h/2),p,scheme.selected(),mode);
      const auto k3=flavour_derivative(flavour_add(y,k2,h/2),p,scheme.selected(),mode);
      const auto k4=flavour_derivative(flavour_add(y,k3,h),p,scheme.selected(),mode);
      y=flavour_add(flavour_add(flavour_add(flavour_add(y,k1,h/6),k2,h/3),k3,h/3),k4,h/6);
    }
    y.b=before_beta+(mode==QEDCouplingMode::fixed?0:p.beta0*p.length);
    flavour_match(y,p,scheme.photon(),mode);
  }
  return y;
}
template<class T> Real flavour_distance(const FlavourState<T>& a,const FlavourState<T>& b,bool relative_coefficients=false) {
  Real d=std::max({std::abs(a.a-b.a),std::abs(a.b-b.b),std::abs(a.integral_a-b.integral_a),std::abs(a.integral_a2-b.integral_a2)});
  for(auto pair:{std::pair<const FlavourColumns<T>*,const FlavourColumns<T>*>{&a.d,&b.d},
                  {&a.e,&b.e},{&a.c1,&b.c1},{&a.c2,&b.c2}})
    for(std::size_t i=0;i<a.d.size();++i)for(unsigned j=0;j<3;++j) {
      Real value=std::abs((*pair.first)[i][j]-(*pair.second)[i][j]);
      if(relative_coefficients && (pair.first==&a.c1 || pair.first==&a.c2))value/=std::max(Real(1),std::abs((*pair.first)[i][j]));
      if(!std::isfinite(value))throw std::runtime_error("Nonfinite flavour evolution");
      d=std::max(d,value);
    }
  return d;
}
template<class T> FlavourState<T> evolve_flavour_columns(const QEDFlavourScheme& s,Real a,
    const FlavourColumns<T>& boundary,const std::vector<FlavourPhase<T>>& phases,
    QEDCouplingMode mode,const MomentEvolutionOptions& options,std::size_t& steps,Real& shift) {
  steps=options.initial_steps;
  auto coarse=integrate_flavour_steps(s,a,boundary,phases,mode,steps);
  while(true) {
    if(steps>options.max_steps/2)throw std::runtime_error("Flavour evolution failed step refinement");
    steps*=2;auto fine=integrate_flavour_steps(s,a,boundary,phases,mode,steps);
    const bool strong=std::any_of(phases.begin(),phases.end(),[](const auto& p){return !p.qcd0.empty();});
    shift=flavour_distance(fine,coarse,strong);if(shift<=options.step_tolerance)return fine;coarse=std::move(fine);
  }
}
QEDMomentEvolution evolve_flavour_moment(Real N,const MomentEvolutionParameters&,const MomentEvolutionOptions&);
struct FlavourXSpaceResult {
  LatticeMatrix weights,propagator,kernel0,kernel1,boundary1,c1,c2;
  QEDXSpaceDiagnostics diagnostics;
  std::vector<std::array<std::vector<Real>,3>> full_weights;
};
FlavourXSpaceResult evolve_flavour_xspace(const MomentEvolutionParameters&,const QEDXSpaceGrid&,const MomentEvolutionOptions&);
}} // namespace leptonff::detail
