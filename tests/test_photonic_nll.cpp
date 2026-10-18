#include "leptonff/photonic_dressing.hpp"
#include "leptonff/nlo_lepton.hpp"
#include "check.hpp"
#include <algorithm>
#include <numeric>
#include <limits>

int main() {
  using namespace leptonff;
  const auto source=nlo_lepton_component(NLOLeptonComponent::photonic);
  const auto conditioned=unitary_photonic_nlo_kernel();
  close(source.moment(1),photonic_source_number_moment(),2e-9L,"source first moment is preserved");
  close(conditioned.moment(1),0,2e-9L,"conditioned source first moment");
  for(Real z:{.01L,.2L,.7L,.99L})
    close(conditioned.resolved_density(z),source.resolved_density(z),0,"conditioning changes only endpoint");
  QEDParameters q;q.mu_frag=80.385L;
  const LLGrid grid{256,40,.08L,1e-10L};
  const PhotonicNLL nll(q,grid);
  close(nll.moment(1),1,2e-12L,"unitary photonic NLL number");
  const auto& u=nll.baseline().evolution();
  const Real a=q.expansion_parameter(),L=2*std::log(q.mu_frag/q.mu_0);
  close(nll.diagnostics().raw_inclusive_factor,1+a*a*L*source.moment(1),1e-12L,"explicit raw inclusive factor");
  close(nll.diagnostics().kernel_number_residual,0,1e-12L,"projected kernel conservation");
  close(nll.diagnostics().evolved_number_residual,0,1e-11L,"signed evolved kernel conservation");
  // Literal circular convolution independently checks the transported insertion.
  for(std::size_t j:{0u,1u,5u,70u,255u}) {
    Real literal=0;
    for(std::size_t k=0;k<grid.bins;++k)
      literal+=u.weights()[k]*nll.kernel_weights()[(j+grid.bins-k)%grid.bins];
    close(nll.evolved_kernel_weights()[j],literal,2e-12L,"transported P1 literal convolution");
  }
  const FourVector mu{50,40,0,30},nu{40,-24,0,-32},other{50,-30,40,0};
  FiducialCuts zc;zc.pt_min=18;zc.mass_min=40;zc.mass_max=85;
  const PhotonicLabDressing ll(q,grid,.01L),matched(q,grid,.01L,PhotonicAccuracy::nll);
  for(bool neutral:{false,true}) {
    const auto before=neutral?ll.fragmentation().z(mu,other,zc):ll.fragmentation().w(mu,nu);
    const auto after=neutral?matched.fragmentation().z(mu,other,zc):matched.fragmentation().w(mu,nu);
    close(after.fixed_order,before.fixed_order,0,"NLL preserves continuum NLO result");
    close(after.lattice_coefficient1,before.lattice_coefficient1,0,"NLL preserves actual first-order overlap");
    Real correction=0,order2=0;
    for(std::size_t i=0;i<grid.bins;++i) {
      const auto pass=[&](Real x,Real y) {
        return neutral?zc.accepts(mu.scaled(x),other.scaled(y)):mu.scaled(x).pt()>25;
      };
      if(pass(u.z(i),1))order2+=L*nll.kernel_weights()[i];
      if(neutral && pass(1,u.z(i)))order2+=L*nll.kernel_weights()[i];
      if(!neutral) {
        if(pass(u.z(i),1))correction+=nll.evolved_kernel_weights()[i];
      } else for(std::size_t j=0;j<grid.bins;++j)if(pass(u.z(i),u.z(j)))
        correction+=nll.evolved_kernel_weights()[i]*u.weights()[j]+u.weights()[i]*nll.evolved_kernel_weights()[j];
    }
    close(after.matched-before.matched,a*a*L*correction,2e-13L,"W/Z exactly one P1 insertion, no boundary cross");
    close(after.coefficient2-before.coefficient2,order2,2e-10L,"W/Z second-order coefficient includes L P1");
  }
  // H1 is already subleading: adding P1 to it would be NNLL.
  close(matched.w_projected_shift(mu,nu),ll.w_projected_shift(mu,nu),0,"hard W retains LL dressing");
  close(matched.z_projected_shift(mu,other,zc),ll.z_projected_shift(mu,other,zc),0,"hard Z retains LL dressing");
  Real previous_w=0,previous_z=0;
  for(Real alpha:{1e-4L,2e-4L}) {
    q.alpha=alpha;const PhotonicLabDressing d(q,grid,.01L,PhotonicAccuracy::nll);
    const auto w=d.fragmentation().w(mu,nu),z=d.fragmentation().z(mu,other,zc);
    const Real aa=q.expansion_parameter();
    const Real rw=w.matched-w.fixed_order-aa*aa*w.coefficient2;
    const Real rz=z.matched-z.fixed_order-aa*aa*z.coefficient2;
    if(previous_w)close(rw/previous_w,8,.12L,"W remainder is cubic");
    if(previous_z)close(rz/previous_z,8,.12L,"Z remainder is cubic");
    previous_w=rw;previous_z=rz;
  }
  q=QEDParameters{};q.alpha=0;
  const PhotonicNLL zero(q,grid);
  close(zero.moment(3),1,1e-15L,"zero coupling identity");
  q=QEDParameters{};q.mu_0=q.mu_frag;
  const PhotonicLabDressing identity(q,grid,.01L,PhotonicAccuracy::nll);
  close(identity.fragmentation().z(mu,other,zc).matched,identity.fragmentation().z(mu,other,zc).fixed_order,
        2e-13L,"zero evolution recovers exact continuum NLO");
  throws<std::invalid_argument>([&]{nll.moment(.5L);},"reject unsupported moment");
  throws<std::invalid_argument>([&]{nll.apply([](Real){return std::numeric_limits<Real>::quiet_NaN();});},"reject nonfinite measurement");
  std::cout<<"Photonic NLL source convention, unitary shape, W/Z matching and re-expansion passed\n";
}
