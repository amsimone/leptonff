#include "leptonff/coupled_lab.hpp"
#include "check.hpp"
#include <algorithm>

int main() {
  using namespace leptonff;
  MomentEvolutionParameters legacy;
  QEDFlavourScheme one;one.fermions={{"mu",-1,legacy.muon_mass,1,0}};
  auto generic=legacy;generic.flavour_scheme=one;
  for(auto mode:{QEDCouplingMode::fixed,QEDCouplingMode::one_loop,QEDCouplingMode::two_loop}) {
    legacy.coupling=generic.coupling=mode;
    const auto old=evolve_qed_moment(3,legacy),now=evolve_qed_moment(3,generic);
    for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j) {
      close(now.fragmentation[i][j],old.fragmentation[i][j],2e-11L,"general solver reduces to legacy FF");
      close(now.evolution[i][j],old.evolution[i][j],2e-11L,"general solver reduces to legacy propagator");
      close(now.coefficient1[i][j],old.coefficient1[i][j],2e-9L,"general C1 reduction");
      close(now.coefficient2[i][j],old.coefficient2[i][j],3e-7L,"general C2 reduction");
    }
    close(run_flavour_coupling(generic,generic.mu_frag).a_final,old.a_final,2e-13L,"analytic generic running vs legacy RK4");
  }
  for(Real factor:{.5L,1.L,2.L}) {
    auto p=generic;p.flavour_scheme=QEDFlavourScheme::leptons(factor);p.coupling=QEDCouplingMode::two_loop;
    p.alpha_msbar_at_mu0=leptonic_alpha_msbar(1/137.035999084L,p.mu_0,*p.flavour_scheme);
    const auto full=evolve_qed_moment(3,p);const auto flow=run_flavour_coupling(p,p.mu_frag);
    close(full.a_final,flow.a_final,2e-13L,"threshold running independent analytic flow");
    close(full.integral_a,flow.integral_a,2e-12L,"threshold integral a");
    close(full.integral_a2,flow.integral_a2,2e-12L,"threshold integral a squared");
    const Real expected=8*std::log(p.mu_frag/p.mu_0)/3+4*std::log(p.mu_frag/1.77693L)/3;
    close(flow.beta0_integral,expected,2e-13L,"coupling threshold logs cancel matching-scale location");
    auto zero=p;zero.alpha_msbar_at_mu0=0;
    close(run_flavour_coupling(zero,p.mu_frag).beta0_integral,expected,2e-13L,"beta integral has a regular zero-coupling limit");
    auto back=p;back.mu_0=p.mu_frag;back.alpha_msbar_at_mu0=flow.a_final*2*std::acos(-1.L);
    close(run_flavour_coupling(back,p.mu_0).a_final,p.alpha_msbar_at_mu0/(2*std::acos(-1.L)),2e-13L,"coupling matching and decoupling inverse");
    const Real threshold=p.flavour_scheme->fermions[2].threshold();auto at=p;at.mu_frag=threshold;
    const auto above=evolve_qed_moment(3,at);
    auto low=at;low.flavour_scheme->fermions.pop_back();const auto below=evolve_qed_moment(3,low);
    const Real a=below.a_final;
    const Real k=flavour_threshold_coefficient(*p.flavour_scheme,2,4,p.flavour_scheme->photon(),threshold).moment(3);
    for(unsigned j=0;j<3;++j) {
      close(above.flavour_fragmentation[4][j],a*k*below.fragmentation[2][j],2e-12L,"new heavy fragmentation finite matching from BELOW photon");
      close(above.flavour_fragmentation[5][j],above.flavour_fragmentation[4][j],1e-14L,"heavy pair charge symmetry");
      close(above.a_final*above.fragmentation[2][j],a*below.fragmentation[2][j],2e-13L,"alpha times photon FF continuous at one-loop matching map");
    }
    Real previous=0;
    for(Real alpha:{.001L,.0005L}) {
      auto small=p;small.alpha_msbar_at_mu0=alpha;const auto r=evolve_qed_moment(3,small);const Real a0=alpha/(2*std::acos(-1.L));
      Real error=0;
      for(std::size_t i=0;i<r.flavour_fragmentation.size();++i)for(unsigned j=0;j<3;++j) {
        const auto selected=p.flavour_scheme->selected();const Real identity=i==selected[j]?1:0;
        error=std::max(error,std::abs(r.flavour_fragmentation[i][j]-identity-a0*r.flavour_coefficient1[i][j]-a0*a0*r.flavour_coefficient2[i][j]));
      }
      if(previous && (previous/error<7 || previous/error>9))throw std::runtime_error("Full matched multi-flavour remainder is not cubic");
      previous=error;
    }
  }
  // FFT kernel projection and all-parent output moments versus the independent
  // continuum moment ODE, plus continuum NLO restoration in the lab operator.
  auto p=generic;p.flavour_scheme=QEDFlavourScheme::leptons();p.coupling=QEDCouplingMode::two_loop;
  // Matching-scale changes can change uncomputed alpha^2 constants, but may
  // not leave an alpha^2 log(Q) term. All full parent rows are checked.
  auto shifted=p;shifted.flavour_scheme=QEDFlavourScheme::leptons(2);
  auto high=p;high.mu_frag*=2;auto shifted_high=shifted;shifted_high.mu_frag*=2;
  const auto c=evolve_qed_moment(3,p),v=evolve_qed_moment(3,shifted);
  const auto hc=evolve_qed_moment(3,high),hv=evolve_qed_moment(3,shifted_high);
  for(std::size_t i=0;i<c.flavour_coefficient2.size();++i)for(unsigned j=0;j<3;++j)
    close(hv.flavour_coefficient2[i][j]-hc.flavour_coefficient2[i][j],
          v.flavour_coefficient2[i][j]-c.flavour_coefficient2[i][j],2e-7L,"no retained alpha2 hard logarithm in threshold-scale shift");
  QEDXSpaceGrid grid;grid.bins=256;
  const QEDXSpaceEvolution x(p,grid);const auto m=evolve_qed_moment(2,p);const auto xm=x.flavour_moment(2);
  for(std::size_t i=0;i<xm.size();++i)for(unsigned j=0;j<3;++j)
    close(xm[i][j],m.flavour_fragmentation[i][j],3e-10L,"full all-parent N2 lattice/continuum agreement");
  const CoupledLabDressing d(SpinDecay::wminus,p,grid);
  const FourVector mu{50,40,0,30},nu{40,-24,0,-32};
  auto oldp=p;oldp.flavour_scheme.reset();const CoupledLabDressing old(SpinDecay::wminus,oldp,grid);
  close(d.fragmentation().w(mu,nu).continuum_coefficient1,old.fragmentation().w(mu,nu).continuum_coefficient1,1e-12L,"muon continuum C1 unchanged by other families");
  close(d.w_projected(mu,nu).beta0_log_hard,run_flavour_coupling(p,invariant_mass(mu,nu)).beta0_integral,1e-13L,"hard matching uses threshold beta integral");
  if(std::abs(d.a_hard(80)-old.a_hard(80))<1e-7L)throw std::runtime_error("Extra charged flavours did not enter hard coupling");
  // Check the ACTUAL lab match, including all hard parents, against its C2
  // in both supported signed fraction measures and all W/Z decay operators.
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector})
    for(auto measure:{FractionMeasure::atoms,FractionMeasure::conservative_cells}) {
      Real previous=0;
      for(Real alpha:{.001L,.0005L}) {
        auto small=p;small.alpha_msbar_at_mu0=alpha;
        const CoupledLabDressing op(kind,small,grid,.01L,measure);
        const bool neutral=kind==SpinDecay::neutral_vector;
        const FourVector partner=neutral?FourVector{50,-30,40,0}:nu;
        const auto ff=neutral?op.fragmentation().z(mu,partner):op.fragmentation().w(mu,partner);
        const auto projection=neutral?op.z_projected(mu,partner):op.w_projected(mu,partner);
        const Real M=invariant_mass(mu,partner);const FourVector parent{mu.energy+partner.energy,mu.px+partner.px,mu.py+partner.py,mu.pz+partner.pz};
        const ParentMeasurements measurement=[&](const FourVector& a,const FourVector& b,const FourVector& g) {
          const Real born=neutral?(FiducialCuts{}.accepts(a,b)?1:0):(a.pt()>25 && std::abs(a.eta())<2.4L?1:0);
          return std::array<Real,3>{born,neutral?op.z_real(a,b,g):op.w_real(a,b,g),
            neutral?op.z_real(a,b,g,{},true):op.w_real(a,b,g,{},true)};
        };
        const auto tensors=decay_hard_geometry_channels(kind,M,.7L,.2L,decay_rotation(.3L,.7L,.9L),parent,measurement);
        SpinTensor density{};for(int i=0;i<3;++i)density[i][i]=1.L/3;
        std::array<Real,3> g;for(int i=0;i<3;++i)g[i]=spin_contract(density,tensors[i]);
        const auto match=CoupledLabDressing::match(ff,projection,g);
        const Real a=op.fragmentation().a_initial();
        const Real error=std::abs(match.matched-match.fixed_order-a*a*match.coefficient2);
        if(previous && (previous/error<7 || previous/error>9))throw std::runtime_error("Threshold-matched lab remainder is not cubic");
        previous=error;
      }
    }
  std::cout<<"Variable-flavour evolution, finite threshold maps, C1/C2 and lab hard coupling passed\n";
}
