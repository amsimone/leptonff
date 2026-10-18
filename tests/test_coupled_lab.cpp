#include "leptonff/coupled_lab.hpp"
#include "leptonff/photon_pair.hpp"
#include "leptonff/vector_current.hpp"
#include "check.hpp"
#include <algorithm>

int main() {
  using namespace leptonff;
  const Real M=91.1876L;
  // An independent photon-pair mass projection uses y=2*k.p_selected/M^2,
  // while this operator uses x=2*E_gamma/M. Their finite remainders differ.
  // Compare the COMPLETE photon channel, for three mass moments, to that
  // separately derived coefficient; do not compare the bare projections.
  for(int power:{1,2,3}) {
    const Real alternative=integrate([&](Real x){return projected_decay_photon_kernel(SpinDecay::neutral_vector,x,M,M)*std::pow(x,power);},0,1).value;
    const Real geometry=integrate([&](Real z){
      return integrate([&](Real c){
        const Real x=1-z,y=x*(1-c)/2,t=x*(1+c)/2;
        const Real real=((1-y)*(1-y)+(1-t)*(1-t))/(y*t)*x/2;
        return real*(std::pow(y,power)+std::pow(t,power)-std::pow(x,power));
      },-1,1).value;
    },0,1).value;
    close(2*alternative+geometry,2*vector_photon_pair_mass_hard_coefficient(M,M).moment(power+1),2e-8L,
          "complete photon channel versus independent pair-mass coefficient");
  }
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector}) {
    for(Real z:{.2L,.6L,.9L}) {
      QEDParameters q;q.mu_frag=M;
      const Real photon_D1=fragmentation_coefficient_order1(Species::muon_minus,Species::photon,q).resolved_density(1-z);
      const Real massive_spectrum=(1+z*z)/(1-z)*(2*std::log(M/q.muon_mass)+std::log(z)-1)
        -(kind==SpinDecay::neutral_vector?0:2*(1+z+z*z)/(3*(1-z)));
      close(photon_D1+projected_decay_photon_kernel(kind,1-z,M,M),massive_spectrum,2e-12L,"photon MSbar subtraction vs massive spectrum");
    }
  }
  // Energy conservation of the complete hard source (unitary convention),
  // and a pointwise photon geometry check independent of the FF prefix code.
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector}) {
    const Real mu=integrate([&](Real z){return projected_decay_hard_kernel(kind,z,M,M)*(z-1);},0,1).value;
    const Real gamma=integrate([&](Real x){return projected_decay_photon_kernel(kind,x,M,M)*x;},0,1).value;
    close(mu+gamma,0,2e-9L,"hard muon plus photon energy sum rule");
  }
  SpinTensor unpolarised{};for(int i=0;i<3;++i)unpolarised[i][i]=1.L/3;
  for(Real z:{.2L,.7L,.95L})for(Real c:{-.6L,.2L,.9L}) {
    const auto f=[&](const FourVector& a,const FourVector& b,const FourVector& k){
      const Real y=std::pow(invariant_mass(a,k)/M,2),t=std::pow(invariant_mass(b,k)/M,2);
      return std::array<Real,3>{0,y*y+t*t,0};
    };
    const auto g=decay_hard_geometry_channels(SpinDecay::neutral_vector,M,z,c,
      decay_rotation(.2L,.7L,1.1L),{M,0,0,0},f);
    const Real x=1-z,y=x*(1-c)/2,t=x*(1+c)/2;
    const Real real=((1-y)*(1-y)+(1-t)*(1-t))/(y*t)*x/2;
    close(spin_contract(unpolarised,g[1]),real*(y*y+t*t-x*x),2e-11L,
      "photon geometry vs independent Dalitz density and mass moment");
  }
  // W energy moments: include the physical angular remainder and every
  // hard parent. Changing mu_frag may leave alpha^2 constants, but its shift
  // must contain no retained log(M/m). Omitting the photon term must fail.
  for(Real N:{2.L,3.L}) {
    const Real geometry=integrate([&](Real z){
      return integrate([&](Real c){
        const auto point=spin_radiative_point(SpinDecay::wminus,M,0,z,c);
        return spin_trace(point.real)/3*(std::pow(2*point.first.energy/M,N-1)-std::pow(z,N-1));
      },-1,1,{1e-7L,1e-7L,8192}).value;
    },0,1,{1e-7L,1e-7L,8192}).value;
    const auto C2=[&](MomentEvolutionParameters p,bool photon) {
      const auto d=evolve_qed_moment(N,p);
      const Real hm=geometry+Distribution{0,{}, {[&](Real z){return projected_decay_hard_kernel(SpinDecay::wminus,z,M,p.mu_frag);}}}.moment(N);
      const Real hg=Distribution{0,[&](Real x){return projected_decay_photon_kernel(SpinDecay::wminus,x,M,p.mu_frag);},{},1}.moment(N);
      return d.coefficient2[0][0]+hm*d.coefficient1[0][0]
        +4*std::log(M/p.mu_0)/3*hm+(photon?hg*d.coefficient1[2][0]:0);
    };
    MomentEvolutionParameters central;central.mu_frag=M;
    auto varied=central;varied.mu_frag=2*M;
    auto light_central=central;light_central.muon_mass/=2;light_central.mu_0=light_central.muon_mass;
    auto light_varied=light_central;light_varied.mu_frag=2*M;
    const Real shift=C2(varied,true)-C2(central,true);
    close(C2(light_varied,true)-C2(light_central,true),shift,2e-7L,
      "W factorisation variation has no retained alpha^2 mass logarithm");
    const Real incomplete=(C2(light_varied,false)-C2(light_central,false))-(C2(varied,false)-C2(central,false));
    if(std::abs(incomplete)<1e-3L)throw std::runtime_error("W missing-photon negative control lost sensitivity");
  }
  MomentEvolutionParameters p;p.mu_frag=M;
  QEDXSpaceGrid grid;grid.bins=256;grid.max_wrap_bound=1e-9L;grid.max_tail_bound=1e-9L;
  const FourVector a{50,40,0,30},b{50,-30,40,0},nu{40,-24,0,-32},photon{40,32,24,0};
  WLabCuts wc;wc.neutrino_pt_min=20;wc.transverse_mass_min=40;
  FiducialCuts zc;zc.pt_min=15;zc.mass_min=20;zc.mass_max=150;
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector}) {
    const CoupledLabDressing d(kind,p,grid,.01L);
    for(auto mode:{QEDCouplingMode::fixed,QEDCouplingMode::one_loop,QEDCouplingMode::two_loop}) {
      auto running=p;running.coupling=mode;const CoupledLabDressing flow(kind,running,grid,.01L);
      for(Real mass:{p.mu_0/2,M/2,M,2*M,13000.L}) {
        running.mu_frag=mass;
        // The moment solver forbids negative evolution intervals; backward
        // running is still defined and covered by the implicit flow above.
        if(mass<p.mu_0)continue;
        close(flow.a_hard(mass),evolve_qed_moment(2,running).a_final,2e-13L,"hard coupling versus independent step-doubled RK4");
      }
    }
    const auto& e=d.evolution();
    // Literal signed fraction sums on all distinct hard-parent pairs, with
    // photon charge symmetry and both charge-flipped primaries retained.
    Real literal=0,first=0;
    for(std::size_t i=0;i<grid.bins && e.z(i)>=.01L;++i) {
      const Real zi=e.z(i),S=e.weights()[0][0][i]/zi,G=e.weights()[2][0][i]/zi;
      if(kind!=SpinDecay::neutral_vector) {
        const auto accept=[&](const FourVector& mu){return mu.pt()>wc.muon_pt_min && std::abs(mu.eta())<wc.muon_eta_max && nu.pt()>wc.neutrino_pt_min && transverse_mass(mu,nu)>wc.transverse_mass_min?1.L:0.L;};
        literal+=S*accept(a.scaled(zi))+G*accept(photon.scaled(zi));
        first+=e.coefficient1()[0][0][i]/zi*accept(a.scaled(zi))+e.coefficient1()[2][0][i]/zi*accept(photon.scaled(zi));
      }else for(std::size_t j=0;j<grid.bins && e.z(j)>=.01L;++j) {
        const Real zj=e.z(j),Sj=e.weights()[0][0][j]/zj;
        literal+=(S*Sj+e.weights()[0][1][i]/zi*e.weights()[0][1][j]/zj)*(zc.accepts(a.scaled(zi),b.scaled(zj))?1:0);
        for(const auto mu:{a,b})literal+=G*(Sj+e.weights()[0][1][j]/zj)*(zc.accepts(photon.scaled(zi),mu.scaled(zj))?1:0);
      }
    }
    close(kind==SpinDecay::neutral_vector?d.z_real(a,b,photon,zc):d.w_real(a,nu,photon,wc),literal,2e-12L,"full real hard measurement vs literal parent sums");
    if(kind!=SpinDecay::neutral_vector)close(d.w_real(a,nu,photon,wc,true),first,2e-12L,"first real coefficient includes photon conversion");
    const auto rotation=decay_rotation(.2L,.7L,1.1L);
    const auto bare=[&](const FourVector& x,const FourVector& y)->Real {
      return kind==SpinDecay::neutral_vector?(zc.accepts(x,y)?1:0):(x.pt()>wc.muon_pt_min && std::abs(x.eta())<wc.muon_eta_max && y.pt()>wc.neutrino_pt_min && transverse_mass(x,y)>wc.transverse_mass_min?1:0);
    };
    const auto f=[&](const FourVector& x,const FourVector& y,const FourVector& k){return std::array<Real,3>{bare(x,y),kind==SpinDecay::neutral_vector?d.z_real(x,y,k,zc):d.w_real(x,y,k,wc),kind==SpinDecay::neutral_vector?d.z_real(x,y,k,zc,true):d.w_real(x,y,k,wc,true)};};
    const auto g=decay_hard_geometry_channels(kind,M,.5L,.3L,rotation,{M,0,0,0},f,1);
    const auto old=decay_hard_geometry(kind,M,.5L,.3L,rotation,{M,0,0,0},bare,1);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)close(std::abs(g[0][i][j]-old[i][j]),0,1e-13L,"bare geometry preserved with shared photon channel");
    for(auto fraction_measure:{FractionMeasure::atoms,FractionMeasure::conservative_cells}) {
    Real oldremainder=0,oldhard=0;
    for(Real alpha:{1e-4L,2e-4L}) {
      auto small=p;small.alpha_msbar_at_mu0=alpha;const CoupledLabDressing low(kind,small,grid,.01L,fraction_measure);
      const auto ff=kind==SpinDecay::neutral_vector?low.fragmentation().z(a,b,zc):low.fragmentation().w(a,nu,wc);
      const auto h=kind==SpinDecay::neutral_vector?low.z_projected(a,b,zc):low.w_projected(a,nu,wc);
      const auto measure=[&](const FourVector& x,const FourVector& y,const FourVector& k){return std::array<Real,3>{bare(x,y),kind==SpinDecay::neutral_vector?low.z_real(x,y,k,zc):low.w_real(x,y,k,wc),kind==SpinDecay::neutral_vector?low.z_real(x,y,k,zc,true):low.w_real(x,y,k,wc,true)};};
      const Real hard_mass=kind==SpinDecay::neutral_vector?invariant_mass(a,b):invariant_mass(a,nu);
      const auto tensor=decay_hard_geometry_channels(kind,hard_mass,.5L,.3L,rotation,{hard_mass,0,0,0},measure,1);
      SpinTensor rho{};rho[0][0]=.7L;rho[1][1]=.2L;rho[2][2]=.1L;
      const auto result=CoupledLabDressing::match(ff,h,{spin_contract(rho,tensor[0]),spin_contract(rho,tensor[1]),spin_contract(rho,tensor[2])});
      const Real a0=alpha/(2*std::acos(-1.L));
      const Real remainder=result.matched-result.fixed_order-a0*a0*result.coefficient2;
      const Real hard=h.a_hard*h.evolved-a0*h.lattice;
      if(oldremainder) {
        close(remainder/oldremainder,8,.15L,"complete coupled C2 leaves cubic remainder, original hard coupling");
        close(hard/oldhard,4,.02L,"hard dressing starts beyond order alpha");
      }
      oldremainder=remainder;oldhard=hard;
    }
    auto zero=p;zero.alpha_msbar_at_mu0=0;const CoupledLabDressing off(kind,zero,grid,.01L,fraction_measure);
    const auto ff=kind==SpinDecay::neutral_vector?off.fragmentation().z(a,b,zc):off.fragmentation().w(a,nu,wc);
    const auto h=kind==SpinDecay::neutral_vector?off.z_projected(a,b,zc):off.w_projected(a,nu,wc);
    close(CoupledLabDressing::match(ff,h,{1,1,1}).matched,ff.born,0,"zero coupling Born identity");
    }
  }
  std::cout<<"Coupled spin W/Z hard parents, independent photon projection, running and C2 checks passed\n";
}
