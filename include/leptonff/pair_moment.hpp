#pragma once
#include "leptonff/moment_evolution.hpp"

namespace leptonff {
struct VectorPairMoment {
  Real N = 0, a_initial = 0, a_hard = 0;
  Real ff = 0, muon_hard = 0, photon_hard = 0, matched = 0;
  Real hard_muon_coefficient = 0, hard_photon_coefficient = 0;
  // matched = 1 + a_initial*C1 + a_initial^2*C2 + O(a_initial^3).
  // C2 omits the unknown constant d2 boundary and hard H2, both beyond NLL.
  Real coefficient1 = 0, coefficient2 = 0;
};
// Angular-inclusive moment of the inclusive opposite-charge pair density,
// weighted by (m_pair^2/M^2)^(N-1), N>=2. Sum pairs from distinct hard parents:
// mu- mu+, gamma mu-, gamma mu+. One observed lepton per selected direction;
// no veto on additional pairs. Unpolarised vector at rest, one lepton family.
// Complete O(alpha) hard channels + coupled NLO evolution/O(alpha) boundary;
// hard a(M) and FF running use the SAME explicit MSbar input at mu_0.
// This supplies a matched NLL moment reference, not a fiducial prediction.
VectorPairMoment vector_pair_moment(Real mass, Real N,
    const MomentEvolutionParameters& parameters = {},
    const MomentEvolutionOptions& options = {});
} // namespace leptonff
