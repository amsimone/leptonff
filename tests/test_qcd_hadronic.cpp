#include "leptonff/hadronic_input.hpp"
#include "leptonff/coupled_lab.hpp"
#include "check.hpp"
#include <algorithm>
int main() {
  using namespace leptonff;
  HadronicAlphaInput input;MomentEvolutionParameters p;
  configure_hadronic_pointlike(p,input,1/137.035999084L);
  throws<std::invalid_argument>([&]{auto invalid=p;invalid.mu_0=3;configure_hadronic_pointlike(invalid,input,1/137.035999084L);},
    "Reject calibration when the input starts above its matching map");
  const auto& scheme=*p.flavour_scheme;
  close(run_flavour_coupling(p,input.reference_scale).a_final*2*std::acos(-1.L),
    hadronic_alpha_msbar_reference(1/137.035999084L,input),2e-15L,"Physical HVP anchor reproduced without double counting");
  for(Real scale:{2.L,5.L,40.L,91.1876L}) {
    const auto active=scheme.active(scale);
    for(std::size_t i=0;i<scheme.size();++i) {
      Real momentum=0;for(std::size_t j=0;j<scheme.size();++j)momentum+=qcd_leading_kernel(scheme,active,i,j).moment(2);
      close(momentum,0,1e-9L,"LO QCD per-parent momentum conservation");
    }
  }
  close(run_qcd_coupling(p,input.reference_scale)*2*std::acos(-1.L),input.alpha_s_reference,2e-15L,"Strong reference input");
  if(run_qcd_coupling(p,2)<=run_qcd_coupling(p,40))throw std::runtime_error("QCD asymptotic freedom failure");
  close(run_qcd_coupling(p,1.99L),0,0,"No perturbative QCD below matching scale");
  auto composed=p;composed.qcd->reference_scale=5;composed.qcd->alpha_s_reference=run_qcd_coupling(p,5)*2*std::acos(-1.L);
  close(run_qcd_coupling(composed,80),run_qcd_coupling(p,80),1e-15L,"Piecewise strong flow composition");
  const auto moment=evolve_qed_moment(3,p);
  auto none=p;none.qcd->alpha_s_reference=0;const auto no_qcd=evolve_qed_moment(3,none);
  if(std::abs(moment.flavour_fragmentation[scheme.gluon()][2])<1e-5L)throw std::runtime_error("Gluon fragmentation channel absent");
  close(no_qcd.flavour_fragmentation[scheme.gluon()][2],0,0,"Zero alpha_s removes gluon channel");
  none.qcd.reset();const auto pure=evolve_qed_moment(3,none);
  for(std::size_t i=0;i<scheme.size();++i)for(unsigned j=0;j<3;++j)
    close(no_qcd.flavour_fragmentation[i][j],pure.flavour_fragmentation[i][j],2e-11L,"Zero alpha_s reduces to pure-QED engine");
  Real previous=0;
  for(Real alpha:{.001L,.0005L}) {
    auto small=p;small.alpha_msbar_at_mu0=alpha;const auto r=evolve_qed_moment(3,small);Real error=0,a=alpha/(2*std::acos(-1.L));
    for(std::size_t i=0;i<scheme.size();++i)for(unsigned j=0;j<3;++j)
      error=std::max(error,std::abs(r.flavour_fragmentation[i][j]-(i==scheme.selected()[j]?1:0)-a*r.flavour_coefficient1[i][j]-a*a*r.flavour_coefficient2[i][j]));
    if(previous && (previous/error<7 || previous/error>9))throw std::runtime_error("QCD-dressed QED expansion remainder is not cubic");
    previous=error;
  }
  QEDXSpaceGrid grid;grid.bins=128;grid.log_span=64;
  const QEDXSpaceEvolution x(p,grid);const auto xm=x.flavour_moment(3);
  for(std::size_t i=0;i<scheme.size();++i)for(unsigned j=0;j<3;++j)
    close(xm[i][j],moment.flavour_fragmentation[i][j],2e-8L,"QCD/QED all-parent lattice moment agreement");
  grid.bins=64;
  for(auto kind:{SpinDecay::wminus,SpinDecay::neutral_vector}) {
    Real old_error=0;
    for(Real alpha:{.001L,.0005L}) {
      auto small=p;small.alpha_msbar_at_mu0=alpha;
      CoupledLabDressing op(kind,small,grid,.01L,FractionMeasure::conservative_cells);
      const FourVector mu{50,40,0,30},other=kind==SpinDecay::wminus?FourVector{40,-24,0,-32}:FourVector{50,-30,40,0};
      const auto ff=kind==SpinDecay::wminus?op.fragmentation().w(mu,other):op.fragmentation().z(mu,other);
      const auto projection=kind==SpinDecay::wminus?op.w_projected(mu,other):op.z_projected(mu,other);
      const Real a=op.fragmentation().a_initial();
      const auto m=CoupledLabDressing::match(ff,projection,{.7L,.7L-.05L*a,-.05L});
      const Real error=std::abs(m.matched-m.fixed_order-a*a*m.coefficient2);
      if(old_error && (old_error/error<7 || old_error/error>9))throw std::runtime_error("Hadronic/QCD lab matching remainder is not cubic");
      old_error=error;
    }
  }
  auto bad=p;bad.qcd->start_scale=.1L;
  throws<std::invalid_argument>([&]{bad.validate();},"Reject nonperturbative QCD start");
  std::cout<<"Hadronic input, quark/gluon evolution, conservation, limits and QED expansions passed\n";
}
