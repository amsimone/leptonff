#pragma once

#include "leptonff/nlo_lepton.hpp"

namespace leptonff {
// Complete one-family, unit-charge, MSbar time-like matrix. Coefficient of
// (alpha/(2*pi))^2, with mu_ren=mu_frag. Inclusive real/virtual pair content.
// Parent -> observed, exactly as for leading_kernel(). No thresholds.
Distribution next_to_leading_kernel(Species parent, Species observed);
enum class NLOPhotonComponent { fermion_to_photon,photon_to_fermion_direct,
  photon_to_fermion_vacuum,photon_to_photon };
// Unit external charge and unit closed-loop charge sum; no multiplicities.
Distribution nlo_photon_component(NLOPhotonComponent);
} // namespace leptonff
