#include "leptonff/pair_moment.hpp"
#include "leptonff/photon_pair.hpp"
#include "leptonff/vector_decay.hpp"
#include "matrix_evolution.hpp"
#include <cmath>
#include <stdexcept>

namespace leptonff {
VectorPairMoment vector_pair_moment(Real M, Real N,
    const MomentEvolutionParameters& p, const MomentEvolutionOptions& options) {
  p.validate(); options.validate();
  if (!std::isfinite(M) || M <= 2*p.muon_mass || !std::isfinite(N) || N < 2)
    throw std::invalid_argument("Pair moment requires M>2m and finite N>=2");
  const auto d = evolve_qed_moment(N,p,options);
  VectorPairMoment r;
  r.N = N; r.a_initial = d.a_initial;
  const Real Lh = 2*(std::log(M)-std::log(p.mu_0));
  std::size_t steps = 0;
  Real shift = 0;
  const MomentMatrix zero{};
  r.a_hard = detail::evolve_generators(zero,zero,d.a_initial,Lh,p.coupling,
                                      options,steps,shift).a;
  r.hard_muon_coefficient = vector_mass_hard_coefficient(M,p.mu_frag)
      .moment(N,options.kernel_integration);
  r.hard_photon_coefficient = vector_photon_pair_mass_hard_coefficient(M,p.mu_frag)
      .moment(N,options.kernel_integration);
  const auto& D = d.fragmentation;
  const Real same = D[0][0], opposite = D[0][1], gamma = D[2][0];
  r.ff = same*same+opposite*opposite;
  r.muon_hard = r.a_hard*r.hard_muon_coefficient*r.ff;
  r.photon_hard = r.a_hard*2*r.hard_photon_coefficient*gamma*(same+opposite);
  r.matched = r.ff+r.muon_hard+r.photon_hard;
  const auto& c1 = d.coefficient1;
  const Real b0 = p.coupling == QEDCouplingMode::fixed ? 0 : 2.0L/3;
  r.coefficient1 = 2*c1[0][0]+r.hard_muon_coefficient;
  r.coefficient2 = 2*d.coefficient2[0][0]+c1[0][0]*c1[0][0]+c1[0][1]*c1[0][1]
      +r.hard_muon_coefficient*(2*c1[0][0]+b0*Lh)
      +2*r.hard_photon_coefficient*c1[2][0];
  return r;
}
} // namespace leptonff
