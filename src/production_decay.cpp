#include "leptonff/production_decay.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <random>
#include <stdexcept>

namespace leptonff {
namespace {
FourVector sum(const FourVector& a,const FourVector& b) {return {a.energy+b.energy,a.px+b.px,a.py+b.py,a.pz+b.pz};}
Real mass(const FourVector& q) {
  for(Real x:{q.energy,q.px,q.py,q.pz})if(!std::isfinite(x))throw std::invalid_argument("Non-finite parent momentum");
  const Real m2=q.energy*q.energy-q.px*q.px-q.py*q.py-q.pz*q.pz;
  if(q.energy<=0 || m2<=0)throw std::invalid_argument("A positive timelike parent is required");
  return std::sqrt(m2);
}
void massless(const FourVector& p) {
  for(Real x:{p.energy,p.px,p.py,p.pz})if(!std::isfinite(x))throw std::invalid_argument("Non-finite Born momentum");
  if(p.energy<=0 || std::abs(p.energy*p.energy-p.px*p.px-p.py*p.py-p.pz*p.pz)>1e-9L*p.energy*p.energy)
    throw std::invalid_argument("Massless Born decay momenta required");
}
void add(SpinTensor& to,const SpinTensor& from,Real weight) {
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)to[i][j]+=weight*from[i][j];
}
std::vector<std::pair<Real,Real>> gauss(unsigned count) {
  std::vector<std::pair<Real,Real>> out(count);const Real pi=std::acos(-1.L);
  for(unsigned i=0;i<(count+1)/2;++i) {
    Real x=std::cos(pi*(i+.75L)/(count+.5L)),derivative=0;
    for(int iteration=0;iteration<100;++iteration) {
      Real p=1,previous=0;
      for(unsigned k=1;k<=count;++k) {const Real next=((2*k-1)*x*p-(k-1)*previous)/k;previous=p;p=next;}
      derivative=count*(x*p-previous)/(x*x-1);
      const Real step=p/derivative;x-=step;
      if(std::abs(step)<2e-16L)break;
      if(iteration==99)throw std::runtime_error("Angular quadrature did not converge");
    }
    const Real weight=2/((1-x*x)*derivative*derivative);
    out[i]={-x,weight};out[count-1-i]={x,weight};
  }
  return out;
}
DecayPrediction from(const std::array<Real,5>& v) {return {v[0],v[1],v[2],v[3],v[4]};}
}
void DecayConfiguration::validate() const {
  if(process!=SpinDecay::wplus && process!=SpinDecay::wminus && process!=SpinDecay::neutral_vector)
    throw std::invalid_argument("Unsupported decay process");
  if(mode!=DecayMode::ff_only && mode!=DecayMode::photonic_ll && mode!=DecayMode::photonic_nll && mode!=DecayMode::coupled_nll)
    throw std::invalid_argument("Unsupported decay mode");
  evolution.validate();
  if(mode==DecayMode::photonic_ll || mode==DecayMode::photonic_nll) {
    LLGrid ll;ll.bins=grid.bins;ll.log_span=grid.log_span;ll.soft_log_cutoff=ll.spacing()/2;ll.validate();
    if(fraction_measure!=FractionMeasure::atoms)throw std::invalid_argument("Photonic mode uses atomic fraction tables");
  }else grid.validate();
  for(Real x:{soft_cutoff,collinear_cutoff})if(!std::isfinite(x) || x<=0 || x>=.01L)throw std::invalid_argument("Invalid finite-geometry cutoff");
  if(!std::isfinite(minimum_z) || minimum_z<=0 || minimum_z>=1 || hard_samples>4096)
    throw std::invalid_argument("Invalid decay sampling/support configuration");
  if((mode==DecayMode::photonic_nll || mode==DecayMode::coupled_nll) && !hard_samples)
    throw std::invalid_argument("Matched NLL requires finite hard sampling");
  if(mode==DecayMode::photonic_nll && !dress_photonic_hard)
    throw std::invalid_argument("Photonic NLL requires boundary and hard dressing");
  if((mode==DecayMode::photonic_ll || mode==DecayMode::photonic_nll) && (evolution.flavour_scheme || evolution.qcd || evolution.hadronic_matching))
    throw std::invalid_argument("Photonic evolution excludes coupled flavour/strong inputs");
}
std::uint64_t radiation_seed(const double* a,const double* b) {
  if(!a || !b)throw std::invalid_argument("Null radiation reference");
  std::uint64_t value=0x7a2f9c4d1103b6e5ULL;
  for(const auto* p:{a,b})for(int k=0;k<4;++k) {std::uint64_t bits;std::memcpy(&bits,p+k,sizeof(bits));value^=bits;value*=1099511628211ULL;}
  return value;
}
FourVector map_massless_parent(const double* p,Real tolerance) {
  if(!p || !std::isfinite(tolerance) || tolerance<0)throw std::invalid_argument("Invalid massless map input");
  const Real norm=std::hypot(std::hypot(Real(p[0]),Real(p[1])),Real(p[2]));
  if(!std::isfinite(norm) || !std::isfinite(p[3]) || p[3]<=0 || std::abs(norm-p[3])>tolerance*p[3])
    throw std::invalid_argument("Host momentum exceeds massless mapping-roundoff tolerance");
  return {norm,p[0],p[1],p[2]};
}
struct DecayCalculator::Impl {
  DecayConfiguration config;
  std::unique_ptr<LabFragmentation> ff;
  std::unique_ptr<PhotonicLabDressing> dressing;
  std::unique_ptr<CoupledLabDressing> coupled;
  explicit Impl(const DecayConfiguration& c):config(c) {
    c.validate();
    const auto& p=c.evolution;
    if(c.mode==DecayMode::photonic_ll || c.mode==DecayMode::photonic_nll) {
      QEDParameters q{p.alpha_msbar_at_mu0,p.muon_mass,p.mu_0,p.mu_frag};
      LLGrid grid;grid.bins=c.grid.bins;grid.log_span=c.grid.log_span;grid.soft_log_cutoff=grid.spacing()/2;
      ff=std::make_unique<LabFragmentation>(LLNonSinglet(q.expansion_parameter()*2*std::log(q.mu_frag/q.mu_0),grid),q,c.minimum_z);
      if(c.dress_photonic_hard && c.hard_samples)
        dressing=std::make_unique<PhotonicLabDressing>(q,grid,c.minimum_z,c.mode==DecayMode::photonic_nll?PhotonicAccuracy::nll:PhotonicAccuracy::ll);
    }else if(c.mode==DecayMode::coupled_nll) {
      coupled=std::make_unique<CoupledLabDressing>(c.process,p,c.grid,c.minimum_z,c.fraction_measure);
      ff=std::make_unique<LabFragmentation>(coupled->fragmentation());
    }else ff=std::make_unique<LabFragmentation>(QEDXSpaceEvolution(p,c.grid),c.minimum_z,c.fraction_measure);
  }
  bool neutral() const {return config.process==SpinDecay::neutral_vector;}
  DecayPrediction prediction(const FourVector& a,const FourVector& b,const LabSelection& s,const std::array<Real,3>& g) const {
    if(ff->a_initial()==0) {const Real born=s.accepts(neutral(),a,b);return {born,born,born,0,0};}
    const auto f=neutral()?ff->z(a,b,s.z):ff->w(a,b,s.w);
    const Real alpha=ff->a_initial(),scale=config.evolution.mu_frag;
    if(coupled) {
      const auto h=neutral()?coupled->z_projected(a,b,s.z):coupled->w_projected(a,b,s.w);
      const auto r=CoupledLabDressing::match(f,h,g);
      return {f.born,r.fixed_order,r.matched,alpha*alpha*r.coefficient2,r.matched-r.fixed_order};
    }
    const Real h=config.hard_samples?(neutral()?z_projected_decay_hard(a,b,scale,s.z):w_projected_decay_hard(a,b,scale,s.w)):0;
    const Real nlo=f.fixed_order+alpha*(h+g[0]);
    if(dressing) {
      const auto boundary=neutral()?dressing->fragmentation().z(a,b,s.z):dressing->fragmentation().w(a,b,s.w);
      const Real shift=neutral()?dressing->z_projected_shift(a,b,s.z):dressing->w_projected_shift(a,b,s.w);
      const Real matched=boundary.matched+alpha*(h+g[0]+shift+g[1]);
      return {f.born,nlo,matched,alpha*alpha*(dressing->nll_evolution()?boundary.coefficient2:f.coefficient2),matched-nlo};
    }
    const Real matched=f.matched+alpha*(h+g[0]);return {f.born,nlo,matched,alpha*alpha*f.coefficient2,matched-nlo};
  }
  void geometry(const FourVector& parent,const std::vector<LabSelection>& selections,std::uint64_t seed,
                const std::function<void(std::size_t,const CoupledGeometry&,Real)>& accept) const {
    if(!config.hard_samples || ff->a_initial()==0)return;
    const Real M=mass(parent),pi=std::acos(-1.L);
    std::vector<ParentMeasurements> measurements;
    for(const auto& selection:selections)measurements.push_back([&,selection](const FourVector& a,const FourVector& b,const FourVector& photon) {
      const Real selected=selection.accepts(neutral(),a,b);
      if(dressing)return std::array<Real,3>{selected,(neutral()?dressing->z_measurement(a,b,selection.z):dressing->w_measurement(a,b,selection.w))-selected,0};
      if(!coupled)return std::array<Real,3>{selected,0,0};
      return std::array<Real,3>{selected,neutral()?coupled->z_real(a,b,photon,selection.z):coupled->w_real(a,b,photon,selection.w),
        neutral()?coupled->z_real(a,b,photon,selection.z,true):coupled->w_real(a,b,photon,selection.w,true)};
    });
    std::mt19937_64 rng(seed);const auto uniform=[&](){return Real((rng()>>11)+.5L)/9007199254740992.L;};
    const Real edge=std::acos(1-config.collinear_cutoff),theta_span=pi-2*edge;
    const Real tmin=std::sqrt(config.soft_cutoff),tspan=std::sqrt(1-1e-8L)-tmin;
    for(unsigned k=0;k<config.hard_samples;++k) {
      const Real t=tmin+tspan*uniform(),z=1-t*t,theta=edge+theta_span*uniform(),c=std::cos(theta);
      const Real orientation=2*uniform()-1,phi=2*pi*uniform(),twist=2*pi*uniform();
      std::vector<std::array<Real,3>> controls;
      for(const auto& s:selections) {
        const Real ref=!neutral() || (M*std::sqrt(z)>s.z.mass_min && M*std::sqrt(z)<s.z.mass_max)?1:0;
        controls.push_back(dressing || !coupled?std::array<Real,3>{ref,0,0}:std::array<Real,3>{ref,ref,0});
      }
      const auto tensors=decay_hard_geometry_linear_batch(config.process,M,z,c,decay_rotation(orientation,phi,twist),parent,measurements,controls);
      const Real jacobian=2*t*tspan*theta_span*std::sin(theta)/config.hard_samples;
      for(std::size_t i=0;i<selections.size();++i)accept(i,tensors[i],jacobian);
    }
  }
};
DecayCalculator::DecayCalculator(const DecayConfiguration& c):impl_(std::make_shared<Impl>(c)) {}
const DecayConfiguration& DecayCalculator::configuration() const {return impl_->config;}
const LabFragmentation& DecayCalculator::fragmentation() const {return *impl_->ff;}
const PhotonicLabDressing* DecayCalculator::photonic() const {return impl_->dressing.get();}
const CoupledLabDressing* DecayCalculator::coupled() const {return impl_->coupled.get();}
std::vector<DecayPrediction> DecayCalculator::cumulative(const BornDecayPoint& point,const std::vector<LabSelection>& selections) const {
  massless(point.first);massless(point.second);
  if(selections.empty())throw std::invalid_argument("Empty decay selection list");
  for(const auto& s:selections){s.w.validate();s.z.validate();}
  const auto q=sum(point.first,point.second);mass(q);
  const auto rest=rest_frame_boost(point.first,q,-1);const Real norm=std::hypot(rest.pt(),rest.pz);
  const auto density=spin_angular_projector(impl_->config.process,{rest.px/norm,rest.py/norm,rest.pz/norm});
  std::vector<std::array<Real,3>> geometry(selections.size());
  impl_->geometry(q,selections,point.radiation_seed,[&](std::size_t i,const CoupledGeometry& g,Real weight) {
    for(int j=0;j<3;++j)geometry[i][j]+=weight*spin_contract(density,g[j]);
  });
  std::vector<DecayPrediction> out;
  for(std::size_t i=0;i<selections.size();++i)out.push_back(impl_->prediction(point.first,point.second,selections[i],geometry[i]));
  return out;
}
DecaySpectrumPrediction DecayCalculator::evaluate(const BornDecayPoint& point,const LabSpectrum& spectrum) const {
  const auto values=cumulative(point,spectrum.selections);DecaySpectrumPrediction out;out.total=values.at(0);
  for(const auto& b:spectrum.bins) {
    std::array<Real,5> bin{};
    for(int j=0;j<5;++j)bin[j]=(b.positive<0?0:values.at(b.positive).values()[j])-(b.negative<0?0:values.at(b.negative).values()[j]);
    out.bins.push_back(from(bin));
  }
  return out;
}
SpinDecayResponse DecayCalculator::response(const FourVector& q,const LabSpectrum& spectrum,std::uint64_t seed,const DecayAngularGrid& grid) const {
  if(grid.cosine_nodes<2 || grid.cosine_nodes>512 || grid.azimuth_nodes<4 || grid.azimuth_nodes>2048)
    throw std::invalid_argument("Invalid decay angular quadrature");
  if(spectrum.selections.empty())throw std::invalid_argument("Empty response spectrum");
  const Real M=mass(q),pi=std::acos(-1.L);std::vector<std::array<SpinTensor,5>> integrals(spectrum.selections.size());
  for(auto [c,w]:gauss(grid.cosine_nodes))for(unsigned k=0;k<grid.azimuth_nodes;++k) {
    const Real phi=2*pi*(k+.5L)/grid.azimuth_nodes,s=std::sqrt(1-c*c);
    const std::array<Real,3> n{s*std::cos(phi),s*std::sin(phi),c};
    const FourVector a{M/2,M*n[0]/2,M*n[1]/2,M*n[2]/2},b{M/2,-a.px,-a.py,-a.pz};
    const auto first=rest_frame_boost(a,q,1),second=rest_frame_boost(b,q,1);
    const auto B=spin_born_tensor(impl_->config.process,n);
    for(std::size_t i=0;i<integrals.size();++i) {
      const auto p=impl_->prediction(first,second,spectrum.selections[i],{}).values();
      for(int j=0;j<5;++j)add(integrals[i][j],B,p[j]*w/(2*grid.azimuth_nodes));
    }
  }
  const Real a=impl_->ff->a_initial(),ah=impl_->coupled?impl_->coupled->a_hard(M):a;
  const auto& p=impl_->config.evolution;
  const Real beta=!impl_->coupled?0:p.flavour_scheme?run_flavour_coupling(p,M).beta0_integral:
    p.coupling==QEDCouplingMode::fixed?0:4*std::log(M/p.mu_0)/3;
  impl_->geometry(q,spectrum.selections,seed,[&](std::size_t i,const CoupledGeometry& g,Real weight) {
    add(integrals[i][1],g[0],weight*a);
    if(impl_->coupled) {add(integrals[i][2],g[1],weight*ah);add(integrals[i][3],g[2],weight*a*a);add(integrals[i][3],g[0],weight*a*a*beta);}
    else {add(integrals[i][2],g[0],weight*a);if(impl_->dressing)add(integrals[i][2],g[1],weight*a);}
  });
  for(auto& p:integrals)for(int i=0;i<3;++i)for(int j=0;j<3;++j)p[4][i][j]=p[2][i][j]-p[1][i][j];
  SpinDecayResponse out;out.total=integrals.at(0);
  for(const auto& b:spectrum.bins) {
    std::array<SpinTensor,5> bin{};
    for(int j=0;j<5;++j) {if(b.positive>=0)add(bin[j],integrals.at(b.positive)[j],1);if(b.negative>=0)add(bin[j],integrals.at(b.negative)[j],-1);}
    out.bins.push_back(bin);
  }
  return out;
}
DecaySpectrumPrediction SpinDecayResponse::contract(const SpinTensor& rho) const {
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(!std::isfinite(rho[i][j].real()) || !std::isfinite(rho[i][j].imag()) || std::abs(rho[i][j]-std::conj(rho[j][i]))>1e-10L)
    throw std::invalid_argument("Production tensor must be finite and Hermitian; signed eigenvalues are allowed");
  const auto project=[&](const std::array<SpinTensor,5>& tensors) {std::array<Real,5> v{};for(int j=0;j<5;++j)v[j]=spin_contract(rho,tensors[j]);return from(v);};
  DecaySpectrumPrediction out;out.total=project(total);for(const auto& b:bins)out.bins.push_back(project(b));return out;
}
} // namespace leptonff
