#include "leptonff/coupled_lab.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
// Energy-weighted linear-in-z projection, including the whole-plus atom.
// Keep one node below positive support, including its left basis. The low-z
// source is never folded onto the retained nodes or normalised to a number.
std::vector<Real> project(const Distribution& d,std::size_t n,Real h) {
  std::vector<Real> w(n,0);const Real z1=std::exp(-h);
  w[0]=d.apply([=](Real z){return z*(z-z1)/(1-z1);},z1);
  w[1]=d.apply([=](Real z){return z*(1-z)/(1-z1);},z1);
  for(std::size_t j=1;j<n;++j) {
    const Real za=std::exp(-h*j),zb=std::exp(-h*(j+1));
    const auto density=[&](Real u){const Real z=std::exp(-u);return z*z*d.resolved_density(z);};
    const Real total=integrate(density,j*h,(j+1)*h).value;
    const Real right=integrate([&](Real u){return density(u)*(za-std::exp(-u))/(za-zb);},j*h,(j+1)*h).value;
    w[j]+=total-right;if(j+1<n)w[j+1]+=right;
  }
  return w;
}
std::vector<Real> transport(std::vector<Real> a,std::vector<Real> b) {
  const auto n=a.size();std::size_t size=1;while(size<2*n)size*=2;
  a.resize(size,0);b.resize(size,0);
  auto out=circular_lattice_convolution(a,b).weights;out.resize(n);return out;
}
void massless(const FourVector& p) {
  const Real square=p.energy*p.energy-p.px*p.px-p.py*p.py-p.pz*p.pz;
  if(!std::isfinite(square) || p.energy<=0 || std::abs(square)>1e-9L*p.energy*p.energy)
    throw std::invalid_argument("Coupled hard measurement needs massless finite parents");
}
Real w_action(const FourVector& mu,const FourVector& nu,const WLabCuts& c,const SignedFractionTable& t) {
  massless(mu);massless(nu);c.validate();
  if(mu.pt()==0 || std::abs(mu.eta())>=c.muon_eta_max || (c.neutrino_pt_min>0 && nu.pt()<=c.neutrino_pt_min))return 0;
  Real lo=c.muon_pt_min/mu.pt();
  if(c.transverse_mass_min>0) {
    const Real mt=transverse_mass(mu,nu);if(mt==0)return 0;
    lo=std::max(lo,std::pow(c.transverse_mass_min/mt,2));
  }
  return t.interval(lo);
}
Real z_action(const FourVector& a,const FourVector& b,const FiducialCuts& c,
              const SignedFractionTable& first,const SignedFractionTable& second,FractionPairWorkspace* workspace=nullptr) {
  massless(a);massless(b);c.validate();
  if(c.pt_min<=0 || c.mass_min<=0)throw std::invalid_argument("Coupled Z hard action needs positive pT and pair mass support");
  if(a.pt()==0 || b.pt()==0 || std::abs(a.eta())>=c.eta_max || std::abs(b.eta())>=c.eta_max)return 0;
  const Real l1=c.pt_min/a.pt(),l2=c.pt_min/b.pt();if(l1>=1 || l2>=1)return 0;
  const Real M=invariant_mass(a,b);if(M==0)return 0;
  const Real lo=std::pow(c.mass_min/M,2),hi=std::pow(c.mass_max/M,2);
  return workspace?workspace->apply(first,second,l1,l2,lo,hi):first.pair(second,l1,l2,lo,hi);
}
}
Real projected_decay_photon_kernel(SpinDecay kind,Real x,Real M,Real scale) {
  return projected_decay_hard_kernel(kind,1-x,M,scale);
}
CoupledLabDressing::CoupledLabDressing(SpinDecay kind,const MomentEvolutionParameters& p,
    const QEDXSpaceGrid& grid,Real minimum,FractionMeasure measure)
    :kind_(kind),evolution_(p,grid),fragmentation_(evolution_,minimum,measure) {
  const Real h=grid.spacing();
  const auto n=std::min(grid.bins,static_cast<std::size_t>(std::ceil(-std::log(minimum)/h))+1);
  const auto trim=[&](const std::vector<Real>& a){return std::vector<Real>(a.begin(),a.begin()+n);};
  std::vector<std::vector<Real>> v(count);
  v[identity].assign(n,0);v[identity][0]=1;
  v[same]=trim(evolution_.weights()[0][0]);v[opposite]=trim(evolution_.weights()[0][1]);
  v[gamma]=trim(evolution_.weights()[2][0]);v[c1]=trim(evolution_.coefficient1()[0][0]);
  v[cg]=trim(evolution_.coefficient1()[2][0]);
  v[h0]=project(Distribution{0,{}, {[=](Real z){return projected_decay_hard_kernel(kind,z,p.mu_frag,p.mu_frag);}}},n,h);
  v[p0]=project(leading_kernel(Species::muon_minus,Species::muon_minus),n,h);
  const auto Hgamma=project(Distribution{0,[=](Real x){return projected_decay_photon_kernel(kind,x,p.mu_frag,p.mu_frag);},{},1},n,h);
  const auto Pgamma=project(leading_kernel(Species::muon_minus,Species::photon),n,h);
  for(const auto& pair:std::array<std::pair<Table,Table>,4>{{{hs,same},{ho,opposite},{hc,c1},{hg,gamma}}})
    v[pair.first]=transport(pair.first==hg?Hgamma:v[h0],v[pair.second]);
  // Photon source contracts gamma -> observed muon; the evolved gamma FF
  // already includes its two-loop feed-back and its one-family boundary.
  v[hcg]=transport(Hgamma,v[cg]);
  for(const auto& pair:std::array<std::pair<Table,Table>,4>{{{ps,same},{po,opposite},{pc,c1},{pg,gamma}}})
    v[pair.first]=transport(pair.first==pg?Pgamma:v[p0],v[pair.second]);
  v[pcg]=transport(Pgamma,v[cg]);
  v[charged]=v[same];for(std::size_t j=0;j<n;++j)v[charged][j]+=v[opposite][j];
  t_.reserve(count);
  for(auto& w:v) {
    for(std::size_t k=0;k<n;++k)w[k]/=std::exp(-h*k);
    t_.emplace_back(h,w,minimum,measure);
  }
}
Real CoupledLabDressing::a_hard(Real M) const {
  if(!std::isfinite(M) || M<=0)throw std::invalid_argument("Invalid original hard mass");
  const auto& p=evolution_.parameters();const Real a0=fragmentation_.a_initial();
  if(p.flavour_scheme)return run_flavour_coupling(p,M).a_final;
  if(a0==0 || p.coupling==QEDCouplingMode::fixed)return a0;
  const Real b=2.L/3,L=2*std::log(M/p.mu_0),u0=1/a0;
  Real u=u0-b*L;
  if(u<=2*std::acos(-1.L))throw std::runtime_error("Hard running coupling outside perturbative range");
  if(p.coupling==QEDCouplingMode::one_loop)return 1/u;
  // Exact implicit two-loop flow in u=1/a, avoiding cancellation of 1/a terms.
  for(int k=0;k<32;++k) {
    const Real f=(u0-u)/b-std::log1p(b*(u0-u)/(b*u+1))/(b*b)-L;
    const Real delta=f*(b*u+1)/u;u+=delta;
    if(!std::isfinite(u) || u<=2*std::acos(-1.L))throw std::runtime_error("Hard running coupling outside perturbative range");
    if(std::abs(delta)<1e-15L*u)return 1/u;
  }
  throw std::runtime_error("Hard two-loop coupling did not converge");
}
CoupledLabProjection CoupledLabDressing::projection(Real M,const std::function<Real(Table,Table)>& f,bool pair) const {
  CoupledLabProjection r;r.a_initial=fragmentation_.a_initial();r.a_hard=a_hard(M);
  const auto& p=evolution_.parameters();const Real log=2*std::log(M/p.mu_frag);
  r.beta0_log_hard=p.flavour_scheme?run_flavour_coupling(p,M).beta0_integral
    :(p.coupling==QEDCouplingMode::fixed?0:4*std::log(M/p.mu_0)/3);
  const auto H=[&](Table h,Table p,Table other){return f(h,other)+log*f(p,other);};
  r.lattice=H(h0,p0,identity);
  r.evolved=H(hs,ps,pair?same:identity)+H(hg,pg,pair?same:identity);
  r.first=H(hc,pc,identity)+H(hcg,pcg,identity);
  if(pair) {
    r.lattice+=f(identity,h0)+log*f(identity,p0);
    r.evolved+=H(ho,po,opposite)+H(hg,pg,opposite)
      +f(same,hs)+f(opposite,ho)+f(same,hg)+f(opposite,hg)
      +log*(f(same,ps)+f(opposite,po)+f(same,pg)+f(opposite,pg));
    r.first+=H(h0,p0,c1)+f(identity,hc)+f(c1,h0)+f(identity,hcg)
      +log*(f(identity,pc)+f(c1,p0)+f(identity,pcg));
  }
  return r;
}
CoupledLabProjection CoupledLabDressing::w_projected(const FourVector& a,const FourVector& b,const WLabCuts& c) const {
  if(kind_==SpinDecay::neutral_vector)throw std::invalid_argument("W action on neutral coupled operator");
  auto r=projection(invariant_mass(a,b),[&](Table i,Table){return w_action(a,b,c,t_[i]);},false);
  r.continuum=w_projected_decay_hard(a,b,evolution_.parameters().mu_frag,c);return r;
}
CoupledLabProjection CoupledLabDressing::z_projected(const FourVector& a,const FourVector& b,const FiducialCuts& c) const {
  if(kind_!=SpinDecay::neutral_vector)throw std::invalid_argument("Z action on charged coupled operator");
  FractionPairWorkspace workspace;
  auto r=projection(invariant_mass(a,b),[&](Table i,Table j){return z_action(a,b,c,t_[i],t_[j],&workspace);},true);
  r.continuum=z_projected_decay_hard(a,b,evolution_.parameters().mu_frag,c);return r;
}
Real CoupledLabDressing::w_real(const FourVector& mu,const FourVector& nu,const FourVector& photon,const WLabCuts& c,bool first) const {
  if(kind_==SpinDecay::neutral_vector)throw std::invalid_argument("W action on neutral coupled operator");
  return w_action(mu,nu,c,t_[first?c1:same])+w_action(photon,nu,c,t_[first?cg:gamma]);
}
Real CoupledLabDressing::z_real(const FourVector& minus,const FourVector& plus,const FourVector& photon,const FiducialCuts& c,bool first) const {
  if(kind_!=SpinDecay::neutral_vector)throw std::invalid_argument("Z action on charged coupled operator");
  FractionPairWorkspace workspace;
  const auto F=[&](const FourVector& a,const FourVector& b,Table i,Table j){return z_action(a,b,c,t_[i],t_[j],&workspace);};
  if(first)return F(minus,plus,c1,identity)+F(minus,plus,identity,c1)
      +F(photon,minus,cg,identity)+F(photon,plus,cg,identity);
  if(t_[same].measure()==FractionMeasure::conservative_cells)
    return F(minus,plus,same,same)+F(minus,plus,opposite,opposite)
      +F(photon,minus,gamma,charged)+F(photon,plus,gamma,charged);
  return F(minus,plus,same,same)+F(minus,plus,opposite,opposite)
      +F(photon,minus,gamma,same)+F(photon,minus,gamma,opposite)
      +F(photon,plus,gamma,same)+F(photon,plus,gamma,opposite);
}
CoupledLabMatch CoupledLabDressing::match(const LabFFResult& ff,const CoupledLabProjection& h,const std::array<Real,3>& g) {
  CoupledLabMatch r;
  r.raw=ff.resummed+h.a_hard*(h.evolved+g[1]);
  r.overlap=h.a_initial*(ff.continuum_coefficient1+h.continuum-ff.lattice_coefficient1-h.lattice);
  r.matched=r.raw+r.overlap;
  r.fixed_order=ff.fixed_order+h.a_initial*(h.continuum+g[0]);
  r.coefficient2=ff.coefficient2+h.first+g[2]+h.beta0_log_hard*(h.lattice+g[0]);
  return r;
}
} // namespace leptonff
