#include "leptonff/moment_evolution.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

namespace leptonff {
namespace {
Real beta(const QEDFlavourScheme& s,Real scale) {
  const auto active=s.active(scale);Real nf=0;
  for(std::size_t f=0;f<active.size();++f)if(active[f] && s.fermions[f].colours==3)++nf;
  return 11.L/2-nf/3;
}
std::vector<Real> edges(const QEDFlavourScheme& s,Real lo,Real hi) {
  std::set<Real> points;
  for(const auto& f:s.fermions)if(f.colours==3 && f.threshold()>std::min(lo,hi) && f.threshold()<std::max(lo,hi))points.insert(f.threshold());
  std::vector<Real> out(points.begin(),points.end());if(hi<lo)std::reverse(out.begin(),out.end());out.push_back(hi);return out;
}
}
Real run_qcd_coupling(const MomentEvolutionParameters& p,Real scale) {
  if(!p.qcd)return 0;
  p.validate();if(!std::isfinite(scale) || scale<=0)throw std::invalid_argument("Invalid QCD scale");
  const auto& q=*p.qcd;if(scale<q.start_scale || q.alpha_s_reference==0)return 0;
  Real b=q.alpha_s_reference/(2*std::acos(-1.L)),previous=q.reference_scale;
  for(Real edge:edges(*p.flavour_scheme,previous,scale)) {
    b/=1+b*beta(*p.flavour_scheme,std::sqrt(edge*previous))*2*std::log(edge/previous);previous=edge;
    if(!std::isfinite(b) || b<0 || b>=1/(2*std::acos(-1.L)))throw std::runtime_error("QCD coupling left perturbative range");
  }
  return b;
}
Real integral_qcd_coupling(const MomentEvolutionParameters& p,Real lower,Real upper) {
  if(!p.qcd)return 0;
  p.validate();if(!std::isfinite(lower) || !std::isfinite(upper) || lower<=0 || upper<lower)throw std::invalid_argument("Invalid QCD integral interval");
  lower=std::max(lower,p.qcd->start_scale);if(upper<=lower)return 0;
  Real b=run_qcd_coupling(p,lower),sum=0,previous=lower;
  for(Real edge:edges(*p.flavour_scheme,lower,upper)) {
    const Real k=beta(*p.flavour_scheme,std::sqrt(edge*previous)),L=2*std::log(edge/previous);
    if(b!=0){const Real denominator=1+b*k*L;sum+=std::log1p(b*k*L)/k;b/=denominator;}previous=edge;
  }
  return sum;
}
} // namespace leptonff
