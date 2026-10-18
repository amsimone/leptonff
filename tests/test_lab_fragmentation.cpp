#include "leptonff/lab_fragmentation.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const SignedFractionTable signed_table(.5L,{.8L,.3L,-.2L,.1L,0,0,0,0},.04L);
  close(signed_table.interval(.2L),.8L+.3L-.2L+.1L,1e-15L,"signed prefix without clipping");
  close(signed_table.interval(.2L,1),.3L-.2L+.1L,1e-15L,"strict endpoint exclusion");
  close(signed_table.interval(std::exp(-.5L)),.8L,1e-15L,"strict lower ownership");
  throws<std::invalid_argument>([&]{signed_table.interval(.01L);},"finite support enforced");

  QEDXSpaceGrid grid;grid.bins=256;
  MomentEvolutionParameters p;p.mu_frag=80.379L;
  const QEDXSpaceEvolution e(p,grid);
  const LabFragmentation ff(e,.01L);
  const FourVector mu{50,40,0,30},nu{40,-24,0,-32};
  WLabCuts wc;wc.muon_pt_min=25;wc.neutrino_pt_min=20;wc.transverse_mass_min=55;
  close(transverse_mass(mu.scaled(.7L),nu),std::sqrt(.7L)*transverse_mass(mu,nu),1e-12L,
      "neutrino unchanged under charged fragmentation");
  const Real lower=std::max(wc.muon_pt_min/mu.pt(),std::pow(wc.transverse_mass_min/transverse_mass(mu,nu),2));
  const auto w=ff.w(mu,nu,wc);
  QEDParameters q{p.alpha_msbar_at_mu0,p.muon_mass,p.mu_0,p.mu_frag};
  const auto d1=fragmentation_coefficient_order1(Species::muon_minus,Species::muon_minus,q);
  close(w.continuum_coefficient1,d1.apply([](Real){return 1;},lower),2e-10L,
      "analytic W threshold action vs independent distribution quadrature");
  Real literal=0,c1=0;
  for(std::size_t k=0;k<grid.bins;++k)if(e.z(k)>lower) {
    literal+=e.weights()[0][0][k]/e.z(k);
    c1+=e.coefficient1()[0][0][k]/e.z(k);
  }
  close(w.resummed,literal,1e-14L,"W original signed lattice action");
  close(w.matched,literal+ff.a_initial()*(w.continuum_coefficient1-c1),1e-14L,
      "W removes actual lattice overlap");
  auto failed=wc;failed.neutrino_pt_min=25;
  close(ff.w(mu,nu,failed).matched,0,0,"neutrino cut independent of FF");
  const FourVector zero_nu_pt{40,0,0,40};
  close(ff.w(mu,zero_nu_pt).born,1,0,"disabled neutrino cut accepts zero pT");

  // Lab frame with recoil; compare the prefix contraction to literal tensor
  // cuts, including charge-flipped pairs and expansion coefficients.
  const FourVector second{50,-30,40,0};
  FiducialCuts zcuts;zcuts.pt_min=18;zcuts.mass_min=40;zcuts.mass_max=85;
  const auto z=ff.z(mu,second,zcuts);
  Real zliteral=0,z2=0;
  for(std::size_t i=0;i<grid.bins && e.z(i)*mu.pt()>zcuts.pt_min;++i)
    for(std::size_t j=0;j<grid.bins && e.z(j)*second.pt()>zcuts.pt_min;++j)
      if(zcuts.accepts(mu.scaled(e.z(i)),second.scaled(e.z(j)))) {
        const Real tilt=e.z(i)*e.z(j);
        zliteral+=(e.weights()[0][0][i]*e.weights()[0][0][j]
                    +e.weights()[0][1][i]*e.weights()[0][1][j])/tilt;
        z2+=(e.coefficient1()[0][0][i]*e.coefficient1()[0][0][j]
            +(i==0?e.coefficient2()[0][0][j]:0)
            +(j==0?e.coefficient2()[0][0][i]:0))/tilt;
      }
  close(z.resummed,zliteral,2e-13L,"Z lab literal same/opposite tensor");
  close(z.coefficient2,z2,1e-10L,"Z second expansion literal tensor");
  const Real m=invariant_mass(mu,second),up=std::pow(zcuts.mass_max/m,2);
  const auto leg=[&](Real pt) {
    const Real lo=std::max(zcuts.pt_min/pt,std::pow(zcuts.mass_min/m,2));
    return d1.apply([&](Real x){return x<up?1.L:0.L;},lo,{},up<1?std::vector<Real>{up}:std::vector<Real>{});
  };
  close(z.continuum_coefficient1,leg(mu.pt())+leg(second.pt()),3e-9L,
      "Z continuum finite mass window vs distribution quadrature");

  p.alpha_msbar_at_mu0=0;
  const LabFragmentation zero(QEDXSpaceEvolution(p,grid),.01L);
  close(zero.w(mu,nu,wc).matched,w.born,1e-14L,"W zero-coupling identity");
  close(zero.z(mu,second,zcuts).matched,z.born,1e-14L,"Z zero-coupling identity");
  Real previous=0;
  for(Real alpha:{.001L,.0005L}) {
    p.alpha_msbar_at_mu0=alpha;
    const LabFragmentation small(QEDXSpaceEvolution(p,grid),.01L);
    const auto r=small.w(mu,nu,wc);
    const Real remainder=std::abs(r.matched-r.fixed_order-small.a_initial()*small.a_initial()*r.coefficient2);
    if(previous && (previous/remainder<7.5L || previous/remainder>8.5L))
      throw std::runtime_error("W alpha-cubed remainder after coefficient2 subtraction");
    previous=remainder;
  }
  // Independent fixed-alpha, photon-only branch has no wrong-sign pairs;
  // removing its actual LL lattice overlap restores the SAME continuum D1.
  q.alpha=1/137.035999084L;q.mu_0=q.muon_mass;q.mu_frag=80.379L;
  LLGrid llgrid;llgrid.bins=512;llgrid.log_span=40;llgrid.soft_log_cutoff=llgrid.spacing()/2;
  const LLNonSinglet ll(q.expansion_parameter()*2*std::log(q.mu_frag/q.mu_0),llgrid);
  const LabFragmentation photonic(ll,q,.01L);
  const auto photonw=photonic.w(mu,nu,wc);
  close(photonw.continuum_coefficient1,w.continuum_coefficient1,1e-12L,"photonic/common continuum FF coefficient");
  close(photonw.matched,photonw.resummed+q.expansion_parameter()*(photonw.continuum_coefficient1-photonw.lattice_coefficient1),
        1e-14L,"photonic actual LL overlap removal");
  close(photonic.opposite().interval(.01L),0,0,"photonic wrong-sign absence");
  q.alpha=0;
  const LabFragmentation zero_photonic(LLNonSinglet(0,llgrid),q,.01L);
  close(zero_photonic.w(mu,nu,wc).matched,w.born,0,"zero-alpha photonic identity");
  throws<std::invalid_argument>([&]{ff.w({50,1,0,0},nu);},"massive parent rejected");
  throws<std::invalid_argument>([&]{auto c=wc;c.muon_pt_min=-1;ff.w(mu,nu,c);},"invalid W cuts");
  std::cout<<"Lab W/Z signed convolution, neutrino kinematics and overlap checks passed\n";
}
