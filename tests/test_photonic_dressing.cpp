#include "leptonff/photonic_dressing.hpp"
#include "leptonff/fragmentation.hpp"
#include "check.hpp"
#include <numeric>
#include <algorithm>

int main() {
  using namespace leptonff;
  QEDParameters q;q.mu_frag=80.385L;
  const LLGrid grid{128,32,.15L,1e-10L};
  const PhotonicLabDressing d(q,grid,.01L);
  const auto& e=d.evolution();const auto& u=e.evolution();
  for(auto k:{std::size_t{1},std::size_t{3},std::size_t{7},std::size_t{50}}) {
    const Real lo=std::exp(-(k+.5L)*grid.spacing());
    const Real hi=std::exp(-std::max(grid.soft_log_cutoff,(k-.5L)*grid.spacing()));
    const Real analytic=integrate([&](Real z){return projected_decay_hard_kernel(SpinDecay::neutral_vector,z,q.mu_frag,q.mu_frag);},lo,hi).value;
    const Real w=integrate([&](Real z){return projected_decay_hard_kernel(SpinDecay::wminus,z,q.mu_frag,q.mu_frag);},lo,hi).value;
    close(d.hard_source()[k],analytic,2e-10L,"hard source vs independent kernel quadrature");
    close(d.hard_source()[k]+d.w_extra_source()[k],w,2e-10L,"W source finite charged-W term");
  }
  close(std::accumulate(d.hard_source().begin(),d.hard_source().end(),Real{0}),0,1e-12L,"unitary hard-source plus atom");
  for(std::size_t j=0;j<grid.bins;++j) {
    Real direct=0;for(std::size_t k=0;k<=j;++k)direct+=d.hard_source()[k]*u.weights()[j-k];
    close(d.evolved_hard()[j],direct,2e-11L,"zero-padded transport vs literal causal convolution");
  }
  const FourVector mu{50,40,0,30},nu{40,-24,0,-32},other{50,-30,40,0};
  WLabCuts wc;wc.neutrino_pt_min=20;wc.transverse_mass_min=55;
  FiducialCuts zc;zc.pt_min=18;zc.mass_min=40;zc.mass_max=85;
  const auto fw=d.fragmentation().w(mu,nu,wc),fz=d.fragmentation().z(mu,other,zc);
  Real literalw=0,literalz=0,literalh=0,literalbareh=0;
  const Real lower=std::max(wc.muon_pt_min/mu.pt(),std::pow(wc.transverse_mass_min/transverse_mass(mu,nu),2));
  for(std::size_t i=0;i<grid.bins;++i) {
    if(u.z(i)>lower)literalw+=u.weights()[i]+q.expansion_parameter()*e.evolved_boundary_weights()[i];
    for(std::size_t j=0;j<grid.bins;++j)if(zc.accepts(mu.scaled(u.z(i)),other.scaled(u.z(j)))) {
      literalz+=u.weights()[i]*u.weights()[j]+q.expansion_parameter()*(e.evolved_boundary_weights()[i]*u.weights()[j]+u.weights()[i]*e.evolved_boundary_weights()[j]);
      literalh+=d.evolved_hard()[i]*u.weights()[j]+u.weights()[i]*d.evolved_hard()[j];
      if(j==0)literalbareh+=d.hard_source()[i];if(i==0)literalbareh+=d.hard_source()[j];
    }
  }
  close(fw.resummed,literalw,2e-13L,"one dressed boundary, W literal lattice");
  close(fz.resummed,literalz,2e-13L,"Z excludes a second boundary insertion");
  // Set hard M=mu_frag by an overall rescaling; isolates the universal h source.
  const Real factor=q.mu_frag/invariant_mass(mu,other);
  const FourVector a{mu.energy*factor,mu.px*factor,mu.py*factor,mu.pz*factor};
  const FourVector b{other.energy*factor,other.px*factor,other.py*factor,other.pz*factor};
  auto scaled=zc;scaled.pt_min*=factor;scaled.mass_min*=factor;scaled.mass_max*=factor;
  close(d.z_projected_shift(a,b,scaled),literalh-literalbareh,2e-10L,"Z H1 on both LL legs, complete lattice overlap removed");
  close(d.w_measurement(mu,nu,wc),u.apply([&](Real z){return z>lower?1:0;}),1e-12L,"geometry uses pure U, not boundary-dressed FF");

  // Independent literal scale term, evaluated at the original hard mass.
  // The off-shell Z parent is ABOVE the observed mass window and migrates in.
  for(bool is_z:{false,true}) {
    const Real factor=is_z?2:1;
    const FourVector first{mu.energy*factor,mu.px*factor,mu.py*factor,mu.pz*factor};
    const auto source_other=is_z?other:nu;
    const FourVector second{source_other.energy*factor,source_other.px*factor,source_other.py*factor,source_other.pz*factor};
    const Real log_scale=2*std::log(invariant_mass(first,second)/q.mu_frag);
    auto source=d.hard_source();
    for(std::size_t k=0;k<grid.bins;++k)
      source[k]+=log_scale*u.generator_weights()[k]+(is_z?0:d.w_extra_source()[k]);
    std::vector<Real> evolved(grid.bins,0);
    for(std::size_t j=0;j<grid.bins;++j)for(std::size_t k=0;k<=j;++k)evolved[j]+=source[k]*u.weights()[j-k];
    Real literal=0;
    if(is_z) {
      for(std::size_t i=0;i<grid.bins;++i)for(std::size_t j=0;j<grid.bins;++j)
        if(zc.accepts(first.scaled(u.z(i)),second.scaled(u.z(j))))
          literal+=evolved[i]*u.weights()[j]+u.weights()[i]*evolved[j]
                   -(j==0?source[i]:0)-(i==0?source[j]:0);
      close(d.z_projected_shift(first,second,zc),literal,2e-10L,"off-shell Z original hard-scale reconstruction");
    } else {
      for(std::size_t i=0;i<grid.bins;++i)if(u.z(i)>lower)literal+=evolved[i]-source[i];
      close(d.w_projected_shift(first,second,wc),literal,2e-10L,"W original hard-scale reconstruction");
    }
  }
  Real oldshift=0,oldremainder=0;
  for(Real alpha:{1e-4L,2e-4L}) {
    q.alpha=alpha;const PhotonicLabDressing small(q,grid,.01L);
    const auto ff=small.fragmentation().z(mu,other,zc);
    const Real remainder=ff.matched-ff.fixed_order-q.expansion_parameter()*q.expansion_parameter()*ff.coefficient2;
    const Real shift=q.expansion_parameter()*small.z_projected_shift(mu,other,zc);
    if(oldshift) {
      close(shift/oldshift,4,.03L,"hard dressing starts at alpha squared");
      close(remainder/oldremainder,8,.1L,"boundary LL second-order extraction");
    }
    oldshift=shift;oldremainder=remainder;
  }
  q=QEDParameters{};q.mu_0=q.mu_frag;
  const PhotonicLabDressing identity(q,grid,.01L);
  close(identity.w_projected_shift(mu,nu,wc),0,0,"W hard identity evolution");
  close(identity.z_projected_shift(mu,other,zc),0,0,"Z hard identity evolution");
  close(identity.fragmentation().z(mu,other,zc).matched,identity.fragmentation().z(mu,other,zc).fixed_order,
        1e-13L,"zero evolution gives exact continuum NLO FF, no d1 squared");
  for(auto kind:{SpinDecay::wminus,SpinDecay::neutral_vector}) {
    const Real M=91.1876L;const auto rotation=decay_rotation(.2L,.7L,1.1L);
    const auto bare=[&](const FourVector& p,const FourVector& r)->Real {return kind==SpinDecay::neutral_vector?(zc.accepts(p,r)?1:0):(p.pt()>25 && std::abs(p.eta())<2.4?1:0);};
    const auto evolved=[&](const FourVector& p,const FourVector& r)->Real {return kind==SpinDecay::neutral_vector?identity.z_measurement(p,r,zc):identity.w_measurement(p,r);};
    const auto pair=decay_hard_geometry_pair(kind,M,.7L,.2L,rotation,{M,0,0,0},bare,evolved,1);
    const auto previous=decay_hard_geometry(kind,M,.7L,.2L,rotation,{M,0,0,0},bare,1);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
      close(std::abs(pair[0][i][j]-previous[i][j]),0,0,"shared-current bare geometry preserved");
      close(std::abs(pair[1][i][j]),0,0,"identity dressed geometry difference vanishes pointwise");
    }
  }
  // Nontrivial evolved geometric action versus a literal sum of separate
  // spin-correlated hard events. This checks BOTH Z fractions and keeps the
  // hard parent mass unchanged when a final-state momentum is fragmented.
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector}) {
    const Real M=91.1876L;const auto rotation=decay_rotation(.2L,.7L,1.1L);
    const auto bare=[&](const FourVector& p,const FourVector& r)->Real {
      return kind==SpinDecay::neutral_vector?(zc.accepts(p,r)?1:0)
        :(p.pt()>25 && std::abs(p.eta())<2.4?1:0);
    };
    const auto evolved=[&](const FourVector& p,const FourVector& r)->Real {
      return kind==SpinDecay::neutral_vector?d.z_measurement(p,r,zc):d.w_measurement(p,r);
    };
    const auto paired=decay_hard_geometry_pair(kind,M,.7L,.2L,rotation,{M,0,0,0},bare,evolved,1);
    const auto base=decay_hard_geometry(kind,M,.7L,.2L,rotation,{M,0,0,0},bare,0);
    SpinTensor literal{};
    for(std::size_t i=0;i<grid.bins;++i) {
      if(u.z(i)*M/2< (kind==SpinDecay::neutral_vector?zc.pt_min:25))continue;
      const std::size_t count=kind==SpinDecay::neutral_vector?grid.bins:1;
      for(std::size_t j=0;j<count;++j) {
        if(kind==SpinDecay::neutral_vector && u.z(j)*M/2<zc.pt_min)continue;
        const Real weight=u.weights()[i]*(kind==SpinDecay::neutral_vector?u.weights()[j]:1);
        const auto fragmented=[&](const FourVector& p,const FourVector& r){
          return bare(p.scaled(u.z(i)),kind==SpinDecay::neutral_vector?r.scaled(u.z(j)):r);
        };
        const auto g=decay_hard_geometry(kind,M,.7L,.2L,rotation,{M,0,0,0},fragmented,0);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)literal[a][b]+=weight*g[a][b];
      }
    }
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
      close(std::abs(paired[1][a][b]-(literal[a][b]-base[a][b])),0,1e-11L,
            "dressed geometric spin tensor versus independent fraction sum");
  }
  q.alpha=0;const PhotonicLabDressing zero(q,grid,.01L);
  close(zero.fragmentation().w(mu,nu,wc).matched,fw.born,0,"zero-coupling dressed W");
  close(zero.fragmentation().z(mu,other,zc).matched,fz.born,0,"zero-coupling dressed Z");
  std::cout<<"Photonic lab boundary/hard dressing, causal transport, identity and expansion checks passed\n";
}
