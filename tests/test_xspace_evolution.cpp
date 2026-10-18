#include "leptonff/ff_decay.hpp"
#include "leptonff/xspace_evolution.hpp"
#include "check.hpp"

namespace {
leptonff::Real distance(const leptonff::MomentMatrix& a,const leptonff::MomentMatrix& b) {
  leptonff::Real result=0;
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j)
    result=std::max(result,std::abs(a[i][j]-b[i][j]));
  return result;
}
}
int main() {
  using namespace leptonff;
  QEDXSpaceGrid grid; grid.bins=128;
  const QEDXSpaceEvolution coarse({},grid);
  grid.bins=256;
  const QEDXSpaceEvolution fine({},grid);
  for(Real N:{2.L,3.L})
    close(distance(fine.moment(N),evolve_qed_moment(N).fragmentation),0,2e-10L,
          "endpoint-interpolated moment recovery");
  const auto reference=evolve_qed_moment(5).fragmentation;
  if(distance(fine.moment(5),reference)>=distance(coarse.moment(5),reference)/2)
    throw std::runtime_error("Higher moment did not improve under grid refinement");
  const auto m2=fine.moment(2);
  for(const auto& row:m2) close(row[0]+row[1]+row[2],1,1e-10L,"coupled x-space momentum");
  const auto& d=fine.diagnostics();
  if(d.weighted_wrap_bound>grid.max_wrap_bound || d.weighted_tail_bound>grid.max_tail_bound)
    throw std::runtime_error("Uncontrolled x-space tail/wrap");
  close(fine.apply(Species::muon_minus,Species::muon_minus,[](Real z){return z;},1),
        fine.weights()[0][0][0],1e-15L,"endpoint atom retained in truncated action");
  auto p=MomentEvolutionParameters{};p.alpha_msbar_at_mu0=0;
  grid.bins=32;
  const QEDXSpaceEvolution zero(p,grid);
  close(zero.apply(Species::muon_minus,Species::muon_minus,[](Real){return 1;},.1L),
        1,1e-14L,"zero-coupling endpoint identity");
  const auto zero_match=decay_xspace_ff_match(91.1876L,{},zero);
  close(zero_match.matched,zero_match.born,1e-14L,"zero-coupling fiducial identity");
  Real previous=0;
  for(Real alpha:{.002L,.001L}) {
    p.alpha_msbar_at_mu0=alpha;
    const QEDXSpaceEvolution small(p,grid);
    const auto match=decay_xspace_ff_match(91.1876L,{},small);
    const Real residual=std::abs(match.matched-match.born
        -alpha/(2*std::acos(-1.0L))*match.continuum_coefficient1);
    if(previous && (previous/residual<3.8L || previous/residual>4.2L))
      throw std::runtime_error("FF fiducial match failed quadratic remainder onset");
    previous=residual;
  }
  throws<std::invalid_argument>([&]{fine.moment(1);},"divergent/unbounded lattice number action");
  throws<std::invalid_argument>([&]{fine.apply(Species::photon,Species::muon_minus,
      [](Real){return 1;},0);},"zero lower bound rejected");
  throws<std::invalid_argument>([]{QEDXSpaceGrid bad;bad.bins=31;QEDXSpaceEvolution e({},bad);},
                                "invalid coupled FFT count");
  throws<std::invalid_argument>([]{QEDXSpaceGrid bad;bad.bins=32;bad.log_span=4;
      QEDXSpaceEvolution e({},bad);},"short wrap domain rejected");
  throws<std::invalid_argument>([]{MomentEvolutionParameters bad;bad.mu_frag=bad.mu_0/2;
      QEDXSpaceEvolution e(bad);},"backward x-space evolution rejected");
  std::cout<<"Coupled x-space endpoint, moments, convergence and FF matching checks passed\n";
}
