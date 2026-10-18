#include "leptonff/photonic_dressing.hpp"
#include "leptonff/vector_decay.hpp"
#include "hard_primitive.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace leptonff {
namespace {
std::vector<Real> add(const std::vector<Real>& a,const std::vector<Real>& b) {
  auto out=a;for(std::size_t k=0;k<out.size();++k)out[k]+=b[k];return out;
}
std::vector<Real> delta(std::size_t n) {std::vector<Real> out(n,0);out[0]=1;return out;}
std::vector<Real> project(const LLGrid& grid,bool extra) {
  grid.validate();std::vector<Real> out(grid.bins,0);
  Real sum=0,correction=0;
  for(std::size_t k=1;k<grid.bins;++k) {
    const Real lo=std::max(grid.soft_log_cutoff,(k-.5L)*grid.spacing()),hi=(k+.5L)*grid.spacing();
    if(lo>=hi)continue;
    out[k]=detail::hard_primitive(lo,extra)-detail::hard_primitive(hi,extra);
    const Real y=out[k]-correction,next=sum+y;correction=(next-sum)-y;sum=next;
  }
  out[0]=-sum;return out;
}
std::vector<Real> transport(const std::vector<Real>& source,const LLNonSinglet& u) {
  if(u.time()==0)return source;
  if(source.size()>(std::size_t{1}<<21))throw std::invalid_argument("Hard causal convolution grid too large");
  // Double-size zero padding: no large-loss cell can wrap to physical support.
  auto a=source,b=u.weights();a.resize(2*a.size(),0);b.resize(a.size(),0);
  auto result=circular_lattice_convolution(a,b).weights;result.resize(source.size());return result;
}
void massless(const FourVector& p) {
  const Real square=p.energy*p.energy-p.px*p.px-p.py*p.py-p.pz*p.pz;
  if(!std::isfinite(square) || p.energy<=0 || std::abs(square)>1e-9L*p.energy*p.energy)
    throw std::invalid_argument("Photonic lab dressing requires finite massless parents");
}
Real w_action(const FourVector& mu,const FourVector& nu,const WLabCuts& cuts,const SignedFractionTable& table) {
  massless(mu);massless(nu);cuts.validate();
  if(mu.pt()==0 || std::abs(mu.eta())>=cuts.muon_eta_max ||
     (cuts.neutrino_pt_min>0 && nu.pt()<=cuts.neutrino_pt_min))return 0;
  Real lo=cuts.muon_pt_min/mu.pt();
  if(cuts.transverse_mass_min>0) {
    const Real mt=transverse_mass(mu,nu);if(mt==0)return 0;
    lo=std::max(lo,std::pow(cuts.transverse_mass_min/mt,2));
  }
  return table.interval(lo);
}
Real z_action(const FourVector& a,const FourVector& b,const FiducialCuts& cuts,
              const SignedFractionTable& first,const SignedFractionTable& second,FractionPairWorkspace* workspace=nullptr) {
  massless(a);massless(b);cuts.validate();
  if(cuts.pt_min<=0)throw std::invalid_argument("Photonic lab Z action requires positive pT support");
  if(a.pt()==0 || b.pt()==0 || std::abs(a.eta())>=cuts.eta_max || std::abs(b.eta())>=cuts.eta_max)return 0;
  const Real l1=cuts.pt_min/a.pt(),l2=cuts.pt_min/b.pt();
  if(l1>=1 || l2>=1)return 0;
  const Real M=invariant_mass(a,b);if(M==0)return 0;
  const Real lo=std::pow(cuts.mass_min/M,2),hi=std::pow(cuts.mass_max/M,2);
  return workspace?workspace->apply(first,second,l1,l2,lo,hi):first.pair(second,l1,l2,lo,hi);
}
}
PhotonicLabDressing::PhotonicLabDressing(const QEDParameters& q,const LLGrid& grid,Real minimum,PhotonicAccuracy accuracy)
    :evolution_(q,grid),
     nll_(accuracy==PhotonicAccuracy::nll?std::optional<PhotonicNLL>(std::in_place,evolution_):std::nullopt),
     fragmentation_(nll_?LabFragmentation(*nll_,minimum):LabFragmentation(evolution_,minimum)),
     hard_(project(grid,false)),w_extra_(project(grid,true)),
     splitting_(evolution_.evolution().generator_weights()),
     evolved_hard_(transport(hard_,evolution_.evolution())),
     u_(grid.spacing(),evolution_.evolution().weights(),minimum),
     delta_(grid.spacing(),delta(grid.bins),minimum),h_(grid.spacing(),hard_,minimum),
     hw_(grid.spacing(),add(hard_,w_extra_),minimum),p_(grid.spacing(),splitting_,minimum),
     hu_(grid.spacing(),evolved_hard_,minimum),
     hwu_(grid.spacing(),transport(add(hard_,w_extra_),evolution_.evolution()),minimum),
     pu_(grid.spacing(),transport(splitting_,evolution_.evolution()),minimum) {}
Real PhotonicLabDressing::w_measurement(const FourVector& a,const FourVector& b,const WLabCuts& c) const {
  return w_action(a,b,c,u_);
}
Real PhotonicLabDressing::z_measurement(const FourVector& a,const FourVector& b,const FiducialCuts& c) const {
  return z_action(a,b,c,u_,u_);
}
Real PhotonicLabDressing::w_projected_shift(const FourVector& a,const FourVector& b,const WLabCuts& c) const {
  if(evolution_.evolution().time()==0)return 0;
  const Real log_scale=2*std::log(invariant_mass(a,b)/evolution_.parameters().mu_frag);
  return w_action(a,b,c,hwu_)-w_action(a,b,c,hw_)
         +log_scale*(w_action(a,b,c,pu_)-w_action(a,b,c,p_));
}
Real PhotonicLabDressing::z_projected_shift(const FourVector& a,const FourVector& b,const FiducialCuts& c) const {
  if(evolution_.evolution().time()==0)return 0;
  const Real log_scale=2*std::log(invariant_mass(a,b)/evolution_.parameters().mu_frag);
  FractionPairWorkspace workspace;
  return z_action(a,b,c,hu_,u_,&workspace)+z_action(a,b,c,u_,hu_,&workspace)
        -z_action(a,b,c,h_,delta_,&workspace)-z_action(a,b,c,delta_,h_,&workspace)
        +log_scale*(z_action(a,b,c,pu_,u_,&workspace)+z_action(a,b,c,u_,pu_,&workspace)
                   -z_action(a,b,c,p_,delta_,&workspace)-z_action(a,b,c,delta_,p_,&workspace));
}
} // namespace leptonff
