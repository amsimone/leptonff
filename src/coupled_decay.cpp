#include "leptonff/coupled_decay.hpp"
#include "leptonff/hard_dressing.hpp"
#include "leptonff/photon_pair.hpp"
#include "leptonff/vector_decay.hpp"
#include "matrix_evolution.hpp"
#include "cut_breakpoints.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace leptonff {
namespace {
std::vector<Real> project(const Distribution& d, std::size_t n, Real h,
                          const IntegrationOptions& options) {
  std::vector<Real> w(n,0);
  const Real z1 = std::exp(-h);
  w[0] = d.apply([=](Real z) { return z*(z-z1)/(1-z1); },z1,options);
  if (n == 1) return w;
  w[1] = d.apply([=](Real z) { return z*(1-z)/(1-z1); },z1,options);
  // Last retained node also receives its left basis from the NEXT cell.
  for (std::size_t j=1;j<n;++j) {
    const Real za=std::exp(-h*j), zb=std::exp(-h*(j+1));
    const auto density = [&](Real u) {
      const Real z=std::exp(-u);
      return z*z*d.resolved_density(z);
    };
    const Real total=integrate(density,j*h,(j+1)*h,options).value;
    const Real right=integrate([&](Real u) {
      const Real z=std::exp(-u);
      return density(u)*(za-z)/(za-zb);
    },j*h,(j+1)*h,options).value;
    w[j]+=total-right;
    if(j+1<n) w[j+1]+=right;
  }
  return w;
}
std::vector<Real> convolve(const std::vector<Real>& a,const std::vector<Real>& b) {
  std::vector<Real> c(a.size(),0);
  for(std::size_t i=0;i<c.size();++i)
    for(std::size_t j=0;j+i<c.size();++j) c[i+j]+=a[i]*b[j];
  return c;
}
Real tensor(const std::vector<Real>& a,const std::vector<std::vector<Real>>& g,
            const std::vector<Real>& b) {
  Real result=0;
  for(std::size_t i=0;i<a.size();++i) for(std::size_t j=0;j<b.size();++j)
    result+=a[i]*g[i][j]*b[j];
  return result;
}
std::vector<Real> sum(const std::vector<Real>& a,const std::vector<Real>& b) {
  auto c=a;
  for(std::size_t i=0;i<c.size();++i) c[i]+=b[i];
  return c;
}
} // namespace

void CoupledDecayOptions::validate() const {
  if(ff_panels<1 || ff_panels>64)
    throw std::invalid_argument("Coupled hard FF panels must be in [1,64]");
  if(ff_fraction_breaks.size()>256) throw std::invalid_argument("Too many FF panel contacts");
  for(Real z:ff_fraction_breaks)
    if(!std::isfinite(z) || z<=0 || z>=1)
      throw std::invalid_argument("FF panel contacts must lie in (0,1)");
  geometry.validate();
  integrate([](Real) { return Real{0}; },0,0,projection);
}

CoupledFiducialDecay::CoupledFiducialDecay(Real M,const FiducialCuts& cuts,
    const QEDXSpaceGrid& grid,const CoupledDecayOptions& options)
    : mass_(M),cuts_(cuts),grid_(grid),options_(options) {
  grid.validate(); options.validate(); cuts.validate();
  MassWindow{cuts.mass_min,cuts.mass_max}.validate(M);
  if(cuts.mass_min<=0) throw std::invalid_argument("Coupled hard matching needs positive mass_min");
  const Real mass_padding=options.average_mass_cells ? std::exp(-grid.spacing()/2) : 1;
  minimum_=std::max(2*cuts.pt_min/M,std::pow(cuts.mass_min/M,2)*mass_padding);
  if(cuts.mass_min>=M || minimum_>=1)
    throw std::invalid_argument("Coupled hard reference needs nonempty kinematic support");
  const Real h=grid.spacing();
  if(minimum_<std::exp(-h*(grid.bins-1)))
    throw std::invalid_argument("Coupled hard support extends below FF grid");
  prefix_=1;
  while(prefix_<grid.bins && std::exp(-h*prefix_)>=minimum_) ++prefix_;
  hm_=project(vector_mass_hard_coefficient(M,M),prefix_,h,options.projection);
  hg_=project(vector_photon_pair_mass_hard_coefficient(M,M),prefix_,h,options.projection);
  pm_=project(leading_kernel(Species::muon_minus,Species::muon_minus),prefix_,h,options.projection);
  pg_=project(leading_kernel(Species::muon_minus,Species::photon),prefix_,h,options.projection);

  fractions_.push_back(1);
  constexpr Real nodes[]{-.861136311594052575L,-.339981043584856265L,
                         .339981043584856265L,.861136311594052575L};
  const auto panels=options.ff_panels;
  const Real loss=std::sqrt(1-minimum_);
  std::vector<Real> interior;
  for(std::size_t k=1;k<panels;++k) interior.push_back(loss*k/panels);
  for(Real z:options.ff_fraction_breaks) if(z>minimum_) interior.push_back(std::sqrt(1-z));
  auto loss_edges=detail::resolved_cut_breakpoints(interior,0,loss);
  loss_edges.insert(loss_edges.begin(),0);loss_edges.push_back(loss);
  for(Real u:loss_edges) panel_edges_.push_back(1-u*u);
  for(std::size_t k=0;k+1<loss_edges.size();++k) {
    const Real lo=panel_edges_[k+1],hi=panel_edges_[k];
    for(Real x:nodes) fractions_.push_back((lo+hi)/2+(hi-lo)*x/2);
  }
  compression_.assign(fractions_.size(),std::vector<Real>(prefix_,0));
  compression_[0][0]=1;
  for(std::size_t j=1;j<prefix_;++j) {
    const Real z=std::exp(-h*j);
    const auto location=std::upper_bound(loss_edges.begin(),loss_edges.end(),std::sqrt(1-z));
    const auto k=std::min(loss_edges.size()-2,static_cast<std::size_t>(location-loss_edges.begin()-1));
    for(std::size_t i=1+4*k;i<1+4*(k+1);++i) {
      Real basis=1;
      for(std::size_t l=1+4*k;l<1+4*(k+1);++l)
        if(l!=i) basis*=(z-fractions_[l])/(fractions_[i]-fractions_[l]);
      compression_[i][j]=basis/z; // undo the energy tilt only on support
    }
  }
  std::vector<Real> active_fractions;
  std::vector<std::vector<Real>> active_compression;
  for(std::size_t i=0;i<fractions_.size();++i) {
    if(std::any_of(compression_[i].begin(),compression_[i].end(),[](Real x){return x!=0;})) {
      active_fractions.push_back(fractions_[i]);
      active_compression.push_back(std::move(compression_[i]));
    }
  }
  fractions_=std::move(active_fractions); compression_=std::move(active_compression);
  build_geometry();
}

void CoupledFiducialDecay::build_geometry() {
  const auto nq=fractions_.size();
  gm_.assign(nq,std::vector<Real>(nq,0)); gg_=gm_;
  for(std::size_t i=0;i<nq;++i) for(std::size_t j=0;j<nq;++j) {
    try {
    if(j>=i) gm_[i][j]=vector_muon_pair_geometry(mass_,cuts_,fractions_[i],fractions_[j],options_.geometry,selection_);
    else gm_[i][j]=gm_[j][i];
    gg_[i][j]=fiducial_vector_photon_pair_hard(mass_,mass_,cuts_,fractions_[i],fractions_[j],options_.geometry,selection_).geometry;
    } catch(const std::runtime_error& error) {
      std::ostringstream context;
      context<<"Hard geometry at fractions "<<fractions_[i]<<','<<fractions_[j]
        <<", selection "<<static_cast<int>(selection_.cut)<<", threshold "<<selection_.threshold<<": "<<error.what();
      throw std::runtime_error(context.str());
    }
  }
  em_.assign(prefix_,0);eg_first_=em_;eg_second_=em_;
  if(options_.exact_endpoint_strips) {
    em_[0]=gm_[0][0];eg_first_[0]=gg_[0][0];eg_second_[0]=gg_[0][0];
    for(std::size_t j=1;j<prefix_;++j) {
    const Real z=std::exp(-grid_.spacing()*j);
    em_[j]=vector_muon_pair_geometry(mass_,cuts_,z,1,options_.geometry,selection_)/z;
    eg_first_[j]=fiducial_vector_photon_pair_hard(mass_,mass_,cuts_,z,1,options_.geometry,selection_).geometry/z;
    eg_second_[j]=fiducial_vector_photon_pair_hard(mass_,mass_,cuts_,1,z,options_.geometry,selection_).geometry/z;
    }
  }
}

CoupledFiducialDecay CoupledFiducialDecay::restricted_mass_window(Real lower,Real upper) const {
  MassWindow{lower,upper}.validate(mass_);
  if(lower<cuts_.mass_min || upper>cuts_.mass_max)
    throw std::invalid_argument("Restricted mass window extends outside the hard operator");
  auto result=*this;
  result.cuts_.mass_min=lower; result.cuts_.mass_max=upper;
  if(lower!=cuts_.mass_min || upper!=cuts_.mass_max) result.build_geometry();
  return result;
}

CoupledFiducialDecay CoupledFiducialDecay::with_lepton_selection(const SingleLeptonSelection& selection) const {
  selection.validate();
  auto result=*this;
  result.selection_=selection;
  if((selection.cut==SingleLeptonCut::pt_above && selection.threshold<=cuts_.pt_min) ||
      (selection.cut==SingleLeptonCut::absolute_eta_below && selection.threshold>=cuts_.eta_max))
    result.selection_={};
  if(result.selection_.cut!=selection_.cut || result.selection_.threshold!=selection_.threshold)
    result.build_geometry();
  return result;
}

std::vector<Real> CoupledFiducialDecay::compress(const std::vector<Real>& weighted) const {
  if(weighted.size()<prefix_) throw std::invalid_argument("Short coupled FF weight array");
  std::vector<Real> c(fractions_.size(),0);
  for(std::size_t i=0;i<c.size();++i) for(std::size_t j=0;j<prefix_;++j) {
    if(!std::isfinite(weighted[j])) throw std::invalid_argument("Non-finite coupled FF weight");
    c[i]+=compression_[i][j]*weighted[j];
  }
  return c;
}

CoupledDecayMatch CoupledFiducialDecay::match(const QEDXSpaceEvolution& e,
    const IntegrationOptions& fixed_options) const {
  if(e.grid().bins!=grid_.bins || e.grid().log_span!=grid_.log_span)
    throw std::invalid_argument("Coupled hard and FF grids differ");
  const auto& p=e.parameters();
  if(mass_<=2*p.muon_mass) throw std::invalid_argument("Vector current below pair threshold");
  CoupledDecayMatch r;
  const Real h=grid_.spacing(), Lf=2*(std::log(mass_)-std::log(p.mu_frag));
  r.a_initial=p.alpha_msbar_at_mu0/(2*std::acos(-1.0L));
  MomentEvolutionOptions running_options;
  running_options.step_tolerance=1e-13L;
  std::size_t steps=0; Real shift=0;
  const Real Lh=2*(std::log(mass_)-std::log(p.mu_0));
  r.a_hard=detail::evolve_generators(MomentMatrix{},MomentMatrix{},r.a_initial,Lh,
      p.coupling,running_options,steps,shift).a;
  const auto truncate=[&](const std::vector<Real>& w) {
    return std::vector<Real>(w.begin(),w.begin()+prefix_);
  };
  const auto S=truncate(e.weights()[0][0]), O=truncate(e.weights()[0][1]),
      G=truncate(e.weights()[2][0]), C=truncate(e.coefficient1()[0][0]),
      CO=truncate(e.coefficient1()[0][1]), CG=truncate(e.coefficient1()[2][0]),
      C2=truncate(e.coefficient2()[0][0]);
  std::vector<Real> I(prefix_,0); I[0]=1;
  const auto qs=compress(S),qo=compress(O),qg=compress(G),qc=compress(C),qcg=compress(CG),qi=compress(I);
  auto hm=hm_,hg=hg_;
  for(std::size_t j=0;j<prefix_;++j) { hm[j]+=2*Lf*pm_[j]; hg[j]+=Lf*pg_[j]; }
  const auto born=[&](const std::vector<Real>& a,const std::vector<Real>& b) {
    auto u=a,v=b;
    for(std::size_t j=0;j<prefix_;++j) { const Real z=std::exp(-h*j); u[j]/=z; v[j]/=z; }
    return decay_lattice_acceptance(mass_,cuts_,h,u,v,options_.average_mass_cells,selection_);
  };
  const auto hp=[&](const std::vector<Real>& H,const std::vector<Real>& a,const std::vector<Real>& b) {
    return born(convolve(H,a),b);
  };
  const auto geometry=[&](const auto& matrix,const auto& first_endpoint,const auto& second_endpoint,
      const auto& a,const auto& b,const auto& qa,const auto& qb) {
    if(!options_.exact_endpoint_strips) return tensor(qa,matrix,qb);
    Real value=a[0]*b[0]*matrix[0][0];
    for(std::size_t j=1;j<prefix_;++j)
      value+=a[j]*b[0]*first_endpoint[j]+a[0]*b[j]*second_endpoint[j];
    for(std::size_t i=1;i<qa.size();++i) for(std::size_t j=1;j<qb.size();++j)
      value+=qa[i]*matrix[i][j]*qb[j];
    return value;
  };
  r.born=born(I,I);
  r.charge_flipped_pairs=born(O,O);
  r.ff=born(S,S)+r.charge_flipped_pairs;
  r.muon_projected=hp(hm,S,S)+hp(hm,O,O);
  r.muon_geometry=geometry(gm_,em_,em_,S,S,qs,qs)+geometry(gm_,em_,em_,O,O,qo,qo);
  r.photon_projected=2*hp(hg,G,sum(S,O));
  r.photon_geometry=2*geometry(gg_,eg_first_,eg_second_,G,sum(S,O),qg,sum(qs,qo));
  const Real H0=hp(hm,I,I)+gm_[0][0];
  r.lattice_coefficient1=2*born(C,I)+H0;
  const Real b0=p.coupling==QEDCouplingMode::fixed ? 0 : 2.0L/3;
  r.coefficient2=2*born(C2,I)+born(C,C)+born(CO,CO)
      +hp(hm,C,I)+hp(hm,I,C)+2*geometry(gm_,em_,em_,C,I,qc,qi)+b0*Lh*H0
      +2*(hp(hg,CG,I)+geometry(gg_,eg_first_,eg_second_,CG,I,qcg,qi));
  const QEDParameters q{p.alpha_msbar_at_mu0,p.muon_mass,p.mu_0,p.mu_frag};
  const auto fixed=fiducial_vector_order1(mass_,cuts_,q,options_.geometry,fixed_options,selection_);
  r.continuum_coefficient1=fixed.full_coefficient;
  r.fixed_order=fixed.rate;
  r.overlap_correction=r.a_initial*(r.continuum_coefficient1-r.lattice_coefficient1);
  r.raw=r.ff+r.a_hard*(r.muon_projected+r.muon_geometry+r.photon_projected+r.photon_geometry);
  r.matched=r.raw+r.overlap_correction;
  for(const auto* w:{&S,&O,&G}) {
    const auto compressed=compress(*w);
    for(Real x:compressed) r.compression_absolute_norm+=std::abs(x);
    for(unsigned k=0;k<4;++k) {
      Real original=0,projected=0;
      for(std::size_t j=0;j<prefix_;++j) original+=(*w)[j]*std::exp(-(Real(k)-1)*h*j);
      for(std::size_t i=0;i<compressed.size();++i) projected+=compressed[i]*std::pow(fractions_[i],k);
      r.compression_moment_residual=std::max(r.compression_moment_residual,std::abs(projected-original));
    }
  }
  for(Real value:{r.matched,r.coefficient2,r.compression_absolute_norm})
    if(!std::isfinite(value)) throw std::runtime_error("Non-finite coupled hard result");
  return r;
}
} // namespace leptonff
