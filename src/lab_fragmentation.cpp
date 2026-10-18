#include "leptonff/lab_fragmentation.hpp"
#include "leptonff/photonic_nll.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <sstream>

namespace leptonff {
namespace {
void momentum(const FourVector& p,bool require_massless=true) {
  for(Real x:{p.energy,p.px,p.py,p.pz})
    if(!std::isfinite(x))throw std::invalid_argument("Non-finite laboratory momentum");
  if(p.energy<=0)throw std::invalid_argument("Laboratory energy must be positive");
  const Real square=p.energy*p.energy-p.px*p.px-p.py*p.py-p.pz*p.pz;
  if(require_massless && std::abs(square)>1e-9L*p.energy*p.energy) {
    std::ostringstream detail;
    detail<<"Laboratory FF action requires massless parent leptons: (E,px,py,pz)=("
      <<p.energy<<','<<p.px<<','<<p.py<<','<<p.pz<<"), m2/E2="<<square/(p.energy*p.energy);
    throw std::invalid_argument(detail.str());
  }
}
std::vector<Real> number(const QEDXSpaceEvolution& e,const LatticeMatrix& m,unsigned i,unsigned j) {
  std::vector<Real> result(m[i][j]);
  for(std::size_t k=0;k<result.size();++k)result[k]/=e.z(k);
  return result;
}
std::vector<Real> ll_order(const LLNonSinglet& e,const QEDParameters& p,int order) {
  const Real L=2*std::log(p.mu_frag/p.mu_0);
  auto result=e.generator_weights();
  if(order==2)result=circular_lattice_convolution(result,result).weights;
  for(auto& x:result)x*=order==1?L:L*L/2;
  return result;
}
std::vector<Real> boundary_order2(const LLBoundaryNonSinglet& e) {
  auto result=ll_order(e.evolution(),e.parameters(),2);
  const auto cross=circular_lattice_convolution(e.evolution().generator_weights(),e.initial_boundary_weights());
  const Real L=2*std::log(e.parameters().mu_frag/e.parameters().mu_0);
  for(std::size_t k=0;k<result.size();++k)result[k]+=L*cross.weights[k];
  return result;
}
bool eta_pass(const FourVector& p,Real cut) {
  return p.pt()>0 && std::abs(p.eta())<cut;
}
}

void WLabCuts::validate() const {
  for(Real x:{muon_pt_min,muon_eta_max,neutrino_pt_min,transverse_mass_min})
    if(!std::isfinite(x))throw std::invalid_argument("Non-finite W cut");
  if(muon_pt_min<=0 || muon_eta_max<=0 || neutrino_pt_min<0 || transverse_mass_min<0)
    throw std::invalid_argument("Invalid W laboratory cuts");
}
Real transverse_mass(const FourVector& mu,const FourVector& nu) {
  momentum(mu,false);momentum(nu,false);
  return std::sqrt(std::max(Real{0},2*(mu.pt()*nu.pt()-mu.px*nu.px-mu.py*nu.py)));
}

SignedFractionTable::SignedFractionTable(Real h,const std::vector<Real>& weights,Real minimum,FractionMeasure measure)
    :minimum_z_(minimum),cell_fraction_(std::tanh(h/2)),measure_(measure) {
  if(!std::isfinite(h) || h<=0 || weights.empty() || !std::isfinite(minimum) || minimum<=0 || minimum>=1)
    throw std::invalid_argument("Invalid signed fraction table");
  if(std::exp(-h*(weights.size()-1))>minimum)
    throw std::invalid_argument("Fraction table does not cover minimum_z");
  prefix_.push_back(0);inverse_prefix_.push_back(0);
  for(std::size_t k=0;k<weights.size();++k) {
    if(!std::isfinite(weights[k]))throw std::invalid_argument("Non-finite signed fraction weight");
    const Real z=std::exp(-h*k);
    const Real lower=k==0?1:z*(1-cell_fraction_),upper=k==0?1:z*(1+cell_fraction_);
    if((measure==FractionMeasure::atoms?z:upper)<minimum)break;
    if(k && weights[k]!=0)resolved_=true;
    lower_.push_back(lower);upper_.push_back(upper);
    z_.push_back(z);weights_.push_back(weights[k]);
    prefix_.push_back(prefix_.back()+weights[k]);
    const Real inverse=k==0 || measure==FractionMeasure::atoms?weights[k]/z
      :weights[k]*std::log1p((upper-lower)/lower)/(upper-lower);
    inverse_prefix_.push_back(inverse_prefix_.back()+inverse);
  }
}
Real SignedFractionTable::interval(Real lower,Real upper) const {
  if(!std::isfinite(lower) || !std::isfinite(upper) || lower<0 || upper<0)
    throw std::invalid_argument("Invalid fraction interval");
  if(lower>=1 || upper<=lower)return 0;
  if(lower<minimum_z_)throw std::invalid_argument("Measurement extends below FF table support");
  if(measure_==FractionMeasure::conservative_cells) {
    Real out=tail(lower,false)-tail(upper,false);
    if(upper==1)out-=weights_[0]; // strict upper excludes the exact atom
    return out;
  }
  const auto first=std::upper_bound(z_.begin(),z_.end(),upper,std::greater<Real>{});
  const auto last=std::lower_bound(z_.begin(),z_.end(),lower,std::greater<Real>{});
  if(first>=last)return 0;
  return prefix_[last-z_.begin()]-prefix_[first-z_.begin()];
}

Real SignedFractionTable::cell_lower(std::size_t i) const { return lower_.at(i); }
Real SignedFractionTable::cell_upper(std::size_t i) const { return upper_.at(i); }
Real SignedFractionTable::tail(Real q,bool inverse) const {
  if(q>=1)return 0;
  if(q<minimum_z_)throw std::invalid_argument("Measurement extends below FF table support");
  const auto& prefix=inverse?inverse_prefix_:prefix_;
  if(measure_==FractionMeasure::atoms) {
    const auto last=std::lower_bound(z_.begin(),z_.end(),q,std::greater<Real>{});
    return prefix[last-z_.begin()];
  }
  // Entire cells above q, followed by at most one partial cell. Prefixes
  // include the atom; harmonic cell edges meet without gaps below cell 1.
  const auto last=std::lower_bound(lower_.begin(),lower_.end(),q,std::greater<Real>{});
  const auto k=static_cast<std::size_t>(last-lower_.begin());
  Real value=prefix[k];
  if(k<z_.size() && upper_[k]>q) {
    const Real density=weights_[k]/(upper_[k]-lower_[k]);
    value+=density*(inverse?std::log1p((upper_[k]-q)/q):upper_[k]-q);
  }
  return value;
}
Real SignedFractionTable::inverse_interval(Real lower,Real upper) const {
  if(!std::isfinite(lower) || !std::isfinite(upper) || lower<0 || upper<0)
    throw std::invalid_argument("Invalid inverse fraction interval");
  if(lower>=1 || upper<=lower)return 0;
  Real value=tail(lower,true)-tail(upper,true);
  if(upper<=1 && measure_==FractionMeasure::atoms) {
    const auto it=std::lower_bound(z_.begin(),z_.end(),upper,std::greater<Real>{});
    if(it!=z_.end() && *it==upper)value-=weights_[it-z_.begin()]/upper;
  }else if(upper==1)value-=weights_[0];
  return value;
}
std::array<Real,2> SignedFractionTable::tails(Real q) const {
  if(q>=1)return {};
  if(q<minimum_z_)throw std::invalid_argument("Measurement extends below FF table support");
  if(measure_==FractionMeasure::atoms) {
    const auto last=std::lower_bound(z_.begin(),z_.end(),q,std::greater<Real>{});
    const auto k=last-z_.begin();return {prefix_[k],inverse_prefix_[k]};
  }
  const auto last=std::lower_bound(lower_.begin(),lower_.end(),q,std::greater<Real>{});
  const auto k=static_cast<std::size_t>(last-lower_.begin());
  std::array<Real,2> out{prefix_[k],inverse_prefix_[k]};
  if(k<z_.size() && upper_[k]>q) {
    const Real density=weights_[k]/(upper_[k]-lower_[k]);
    out[0]+=density*(upper_[k]-q);out[1]+=density*std::log1p((upper_[k]-q)/q);
  }
  return out;
}
Real SignedFractionTable::integrated_tail(Real product,Real a,Real b,Real lower) const {
  if(product>=b)return 0;
  if(product==0)return (b-a)*tail(lower,false);
  const auto L=tails(std::max(lower,product/b)),U=tails(std::max(lower,product/a));
  // Exact integral_a^b C(max(lower,product/z)) dz by reversing the order.
  return (b-a)*U[0]+b*(L[0]-U[0])-product*(L[1]-U[1]);
}
Real SignedFractionTable::pair(const SignedFractionTable& second,Real l1,Real l2,Real lo,Real hi) const {
  for(Real x:{l1,l2,lo,hi})if(!std::isfinite(x) || x<0)throw std::invalid_argument("Invalid fraction-pair bounds");
  if(l1>=1 || l2>=1 || hi<=lo || lo>=1)return 0;
  if(l1<minimum_z_ || l2<second.minimum_z_)throw std::invalid_argument("Pair measurement below table support");
  if(measure_==FractionMeasure::conservative_cells) {
    if(!resolved_)return weights_[0]*second.interval(std::max(l2,lo),hi);
    if(!second.resolved_)return second.weights_[0]*interval(std::max(l1,lo),hi);
    l1=std::max(l1,lo);
  }
  if(measure_==FractionMeasure::atoms) {
    Real out=0;
    for(std::size_t k=0;k<z_.size() && z_[k]>l1;++k)
      out+=weights_[k]*second.interval(std::max(l2,lo/z_[k]),hi/z_[k]);
    return out;
  }
  Real out=weights_[0]*second.interval(std::max(l2,lo),hi);
  for(std::size_t k=1;k<z_.size() && upper_[k]>l1;++k) {
    if(weights_[k]==0)continue;
    const Real a=std::max(l1,lower_[k]),b=upper_[k];
    out+=weights_[k]/(upper_[k]-lower_[k])*
      (second.integrated_tail(lo,a,b,l2)-second.integrated_tail(hi,a,b,l2));
  }
  return out;
}

LabFragmentation::LabFragmentation(const QEDXSpaceEvolution& e,Real minimum,FractionMeasure measure)
    :a_(e.parameters().alpha_msbar_at_mu0/(2*std::acos(-1.L))),
     mass_log_(2*std::log(e.parameters().mu_frag/e.parameters().muon_mass)),
     same_(e.grid().spacing(),number(e,e.weights(),0,0),minimum,measure),
     opposite_(e.grid().spacing(),number(e,e.weights(),0,1),minimum,measure),
     c1_(e.grid().spacing(),number(e,e.coefficient1(),0,0),minimum,measure),
     c2_(e.grid().spacing(),number(e,e.coefficient2(),0,0),minimum,measure) {}

LabFragmentation::LabFragmentation(const LLNonSinglet& e,const QEDParameters& p,Real minimum)
    :a_(p.expansion_parameter()),mass_log_(2*std::log(p.mu_frag/p.muon_mass)),
     same_(e.grid().spacing(),e.weights(),minimum),
     opposite_(e.grid().spacing(),std::vector<Real>(e.weights().size(),0),minimum),
     c1_(e.grid().spacing(),ll_order(e,p,1),minimum),
     c2_(e.grid().spacing(),ll_order(e,p,2),minimum) {
  if(std::abs(e.time()-p.expansion_parameter()*2*std::log(p.mu_frag/p.mu_0))>1e-14L ||
     p.mu_0!=p.muon_mass)
    throw std::invalid_argument("Photonic LL lab action requires consistent time and mu0=m_mu");
}

LabFragmentation::LabFragmentation(const LLBoundaryNonSinglet& e,Real minimum)
    :a_(e.parameters().expansion_parameter()),
     mass_log_(2*std::log(e.parameters().mu_frag/e.parameters().muon_mass)),
     same_(e.evolution().grid().spacing(),e.evolution().weights(),minimum),
     opposite_(e.evolution().grid().spacing(),std::vector<Real>(e.evolution().weights().size(),0),minimum),
     c1_(e.evolution().grid().spacing(),e.first_order_weights(),minimum),
     c2_(e.evolution().grid().spacing(),boundary_order2(e),minimum),
     boundary_(std::in_place,e.evolution().grid().spacing(),e.evolved_boundary_weights(),minimum),
     initial_boundary_(std::in_place,e.evolution().grid().spacing(),e.initial_boundary_weights(),minimum) {}

LabFragmentation::LabFragmentation(const PhotonicNLL& e,Real minimum)
    :a_(e.baseline().parameters().expansion_parameter()),
     mass_log_(2*std::log(e.baseline().parameters().mu_frag/e.baseline().parameters().muon_mass)),
     same_(e.baseline().evolution().grid().spacing(),e.baseline().evolution().weights(),minimum),
     opposite_(e.baseline().evolution().grid().spacing(),std::vector<Real>(e.kernel_weights().size(),0),minimum),
     c1_(e.baseline().evolution().grid().spacing(),e.baseline().first_order_weights(),minimum),
     c2_(e.baseline().evolution().grid().spacing(),e.coefficient2_weights(),minimum),
     boundary_(std::in_place,e.baseline().evolution().grid().spacing(),e.finite_insertion_weights(),minimum),
     initial_boundary_(std::in_place,e.baseline().evolution().grid().spacing(),e.baseline().initial_boundary_weights(),minimum) {}

Real LabFragmentation::continuum_interval(Real lower,Real upper) const {
  if(lower>=1 || upper<=lower)return 0;
  if(lower<same_.minimum_z())throw std::invalid_argument("Continuum measurement below FF support");
  // D1 = [(1+z^2)/(1-z)*(log(mu_frag^2/m_mu^2)-1-2log(1-z))]+.
  // Integral from zero to x, analytic in log(1-x); finite lower-bound term
  // of the whole-plus distribution. The endpoint is included only if upper>1.
  const auto primitive=[&](Real x) {
    const Real t=1-x,l=std::log1p(-x);
    const Real j=1.75L-(l*l-2*(t*l-t)+t*t*(l/2-.25L));
    return (mass_log_-1)*(-2*l-x-x*x/2)-2*j;
  };
  if(upper==1)throw std::invalid_argument("Strict continuum endpoint veto is infrared unsafe");
  if(upper>1)return -primitive(lower);
  return primitive(upper)-primitive(lower);
}
LabFFResult LabFragmentation::finish(LabFFResult r) const {
  r.fixed_order=r.born+a_*r.continuum_coefficient1;
  r.overlap_correction=a_*(r.continuum_coefficient1-r.lattice_coefficient1);
  r.matched=r.resummed+r.overlap_correction;
  return r;
}
LabFFResult LabFragmentation::w(const FourVector& mu,const FourVector& nu,const WLabCuts& cuts) const {
  momentum(mu);momentum(nu);cuts.validate();
  LabFFResult r;
  if(!eta_pass(mu,cuts.muon_eta_max) ||
      (cuts.neutrino_pt_min>0 && nu.pt()<=cuts.neutrino_pt_min))return finish(r);
  Real lower=cuts.muon_pt_min/mu.pt();
  if(cuts.transverse_mass_min>0) {
    const Real mt=transverse_mass(mu,nu);
    if(mt==0)return finish(r);
    lower=std::max(lower,std::pow(cuts.transverse_mass_min/mt,2));
  }
  r.born=lower<1?1:0;
  r.resummed=same_.interval(lower);
  if(boundary_)r.resummed+=a_*boundary_->interval(lower);
  r.lattice_coefficient1=c1_.interval(lower);
  r.coefficient2=c2_.interval(lower);
  r.continuum_coefficient1=continuum_interval(lower);
  return finish(r);
}

LabFFResult LabFragmentation::z(const FourVector& first,const FourVector& second,const FiducialCuts& cuts) const {
  momentum(first);momentum(second);cuts.validate();
  if(cuts.pt_min<=0)throw std::invalid_argument("Coupled lab Z action requires positive pT cut");
  LabFFResult r;
  if(!eta_pass(first,cuts.eta_max) || !eta_pass(second,cuts.eta_max))return finish(r);
  const Real l1=cuts.pt_min/first.pt(),l2=cuts.pt_min/second.pt();
  if(l1>=1 || l2>=1)return finish(r);
  const Real mass=invariant_mass(first,second);
  if(mass==0)return finish(r);
  const Real ml=std::pow(cuts.mass_min/mass,2),mu=std::pow(cuts.mass_max/mass,2);
  r.born=(ml<1 && mu>1)?1:0;
  const auto tensor=[&](const SignedFractionTable& a,const SignedFractionTable& b) {
    return a.pair(b,l1,l2,ml,mu);
  };
  r.resummed=tensor(same_,same_)+tensor(opposite_,opposite_);
  if(boundary_)r.resummed+=a_*(tensor(*boundary_,same_)+tensor(same_,*boundary_));
  const Real b1=std::max(l1,ml),b2=std::max(l2,ml);
  r.lattice_coefficient1=c1_.interval(b1,mu)+c1_.interval(b2,mu);
  r.coefficient2=c2_.interval(b1,mu)+c2_.interval(b2,mu)+tensor(c1_,c1_);
  if(initial_boundary_)r.coefficient2-=tensor(*initial_boundary_,*initial_boundary_);
  r.continuum_coefficient1=continuum_interval(b1,mu)+continuum_interval(b2,mu);
  return finish(r);
}
} // namespace leptonff
