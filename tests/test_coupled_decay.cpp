#include "leptonff/coupled_decay.hpp"
#include "leptonff/pair_moment.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const Real M=91.1876L;
  QEDXSpaceGrid grid; grid.bins=256;
  CoupledDecayOptions options; options.ff_panels=2;
  const CoupledFiducialDecay hard(M,FiducialCuts{},grid,options);
  MomentEvolutionParameters p;
  MomentEvolutionOptions evolution_options;
  evolution_options.step_tolerance=1e-13L;
  Real previous=0;
  for(Real alpha:{.00025L,.0005L}) {
    p.alpha_msbar_at_mu0=alpha;
    const auto e=QEDXSpaceEvolution(p,grid,evolution_options);
    const auto r=hard.match(e);
    const Real a=r.a_initial;
    const Real remainder=r.matched-r.born-a*r.continuum_coefficient1-a*a*r.coefficient2;
    if(previous!=0) close(remainder/previous,8,.15L,"coupled fiducial cubic remainder");
    previous=remainder;
    close(r.fixed_order,r.born+a*r.continuum_coefficient1,1e-13L,"complete continuum Oalpha target");
    close(r.compression_moment_residual,0,2e-12L,"linear FF projection moments through cubic");
  }
  p.alpha_msbar_at_mu0=1/137.035999084L;
  const auto central=hard.match(QEDXSpaceEvolution(p,grid));
  close(central.a_hard,vector_pair_moment(M,2,p).a_hard,1e-14L,"common hard coupling with moment reference");
  for(auto mode:{QEDCouplingMode::fixed,QEDCouplingMode::two_loop}) {
    p.coupling=mode;
    const auto base=hard.match(QEDXSpaceEvolution(p,grid));
    p.mu_frag=2*M;
    const auto varied=hard.match(QEDXSpaceEvolution(p,grid));
    close(varied.continuum_coefficient1,base.continuum_coefficient1,2e-7L,"full fiducial Oalpha scale cancellation");
    const Real shift=varied.coefficient2-base.coefficient2;
    p.muon_mass/=2; p.mu_0=p.muon_mass;
    const auto small_varied=hard.match(QEDXSpaceEvolution(p,grid));
    p.mu_frag=M;
    const auto small_base=hard.match(QEDXSpaceEvolution(p,grid));
    close(small_varied.coefficient2-small_base.coefficient2,shift,2e-6L,
          "fiducial scale variation has no alpha2 mass logarithm");
    p.muon_mass*=2; p.mu_0=p.muon_mass;
  }
  p.alpha_msbar_at_mu0=0;
  close(hard.match(QEDXSpaceEvolution(p,grid)).matched,central.born,1e-14L,"coupled zero-alpha identity");
  auto other=grid; other.bins*=2;
  throws<std::invalid_argument>([&]{ hard.match(QEDXSpaceEvolution(p,other)); },"mismatched hard/evolution grid");
  throws<std::invalid_argument>([&]{ CoupledFiducialDecay(M,{25,2.4L,0,120},grid); },"unresolved pair mass");
  std::cout<<"Coupled fiducial hard channels, running and expansion checks passed\n";
}
