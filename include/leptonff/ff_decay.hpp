#pragma once
#include "leptonff/xspace_evolution.hpp"
#include "leptonff/decay.hpp"

namespace leptonff {
struct XSpaceFFDecayMatch {
  Real born = 0, resummed = 0, charge_flipped_pairs = 0;
  Real lattice_coefficient1 = 0, continuum_coefficient1 = 0;
  Real additive_correction = 0, matched = 0;
};
// Born two-direction, inclusive identified opposite-charge pair density:
// D_same x D_same + D_opposite x D_opposite, one observed lepton per direction.
// Exactly O(alpha) FF matched, with no hard H1 and no exclusive extra-pair veto.
// This is FF-layer validation, NOT a full fiducial NLL cross section.
XSpaceFFDecayMatch decay_xspace_ff_match(Real mass,const FiducialCuts& cuts,
    const QEDXSpaceEvolution& evolution,const IntegrationOptions& options = {});
} // namespace leptonff
