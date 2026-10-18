#include "leptonff/decay.hpp"
#include "leptonff/evolution.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  using S = Species;
  const LLNonSinglet small(.001L);
  const FiducialCuts fiducial;
  const Real born = isotropic_acceptance(91.1876L,1,1,fiducial);
  const Real fiducial_coefficient = 2*leading_kernel(S::muon_minus,S::muon_minus).apply(
      [&](Real z) { return isotropic_acceptance(91.1876L,z,1,fiducial); },
      0,{},one_leg_cut_breakpoints(91.1876L,fiducial));
  close((decay_ll_acceptance(91.1876L,fiducial,small)-born)/.001L,
        fiducial_coefficient,.01L,"LL fiducial first-order re-expansion");

  // Independent literal tensor sum verifies the fast min(z1,z2)/mass-window
  // reduction, including parents that migrate in from above the mass window.
  const LLGrid small_grid{512, 32, .04L, 1e-10L};
  const LLNonSinglet fixture(.015L, small_grid);
  auto eta_limited = fiducial;
  eta_limited.eta_max = .3L;
  auto mass_only = fiducial;
  mass_only.pt_min = 0;
  for (Real mass : {91.1876L, 130.0L})
    for (const auto& cuts : {fiducial,eta_limited,mass_only}) {
    Real direct = 0;
    for (std::size_t i = 0; i < fixture.weights().size(); ++i)
      for (std::size_t j = 0; j < fixture.weights().size(); ++j)
        direct += fixture.weights()[i]*fixture.weights()[j]
            *isotropic_acceptance(mass,fixture.z(i),fixture.z(j),cuts);
    close(decay_ll_acceptance(mass,cuts,fixture), direct, 1e-11L,
          "prefix two-leg convolution vs literal tensor cut");
    }
  FiducialCuts inclusive;
  inclusive.pt_min = 0;
  inclusive.mass_min = 0;
  inclusive.eta_max = 20;
  close(decay_ll_acceptance(91.1876L,inclusive,fixture), 1, 1e-11L,
        "LL inclusive two-leg lepton number");

  std::cout << "LL moments, composition, expansion, endpoints and decay checks passed\n";
}
