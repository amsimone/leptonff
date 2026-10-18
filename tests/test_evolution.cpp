#include "leptonff/decay.hpp"
#include "leptonff/evolution.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  using S = Species;
  const LLGrid coarse{4096, 32, .005L, 1e-10L};
  const LLNonSinglet identity(0, coarse);
  close(identity.apply([](Real z) { return z*z; }), 1, 1e-15L,
        "equal-scale LL identity");
  for (Real N : {.5L, 1.0L, 1.5L, 2.0L, 3.0L, 10.0L})
    close(non_singlet_kernel_moment(N),
          leading_kernel(S::muon_minus,S::muon_minus).moment(N), 2e-9L,
          "analytic Mellin kernel vs independent plus action");
  close(non_singlet_kernel_moment(2), -4.0L/3, 1e-12L, "LL P0 N2");
  close(non_singlet_kernel_moment(3), -25.0L/12, 1e-12L, "LL P0 N3");

  constexpr Real time = .015L;
  const LLNonSinglet evolved(time);
  const LLGrid finer{524288, 32, 1e-4L, 1e-10L};
  const LLNonSinglet refined(time, finer);
  close(evolved.moment(1), 1, 1e-12L, "LL lepton number conservation");
  close(evolved.weights()[0], evolved.diagnostics().no_resolved_emission_weight,
        1e-11L, "explicit no-resolved-emission endpoint atom");
  if (evolved.diagnostics().wrap_probability_bound > 1e-10L ||
      evolved.diagnostics().negative_weight_mass > 1e-10L)
    throw std::runtime_error("LL alias/FFT diagnostics failed");
  for (Real N : {2.0L, 3.0L, 10.0L}) {
    const Real exact = evolved.exact_moment(N);
    const Real error = std::abs(evolved.moment(N)-exact);
    const Real refined_error = std::abs(refined.moment(N)-exact);
    close(evolved.moment(N), exact,
          (N-1)*evolved.diagnostics().mean_log_loss_discretisation_bound+1e-10L,
          "LL numerical moments within discretisation bound");
    if (!(refined_error < error))
      throw std::runtime_error("LL moment failed cutoff/grid convergence");
  }

  // A genuine multi-step convolution check in log-loss space, including the
  // endpoint atom, rather than just checking the exponential formula twice.
  const LLNonSinglet first(.007L, coarse), second(.011L, coarse), joined(.018L, coarse);
  for (std::size_t j = 0; j < 64; ++j) {
    Real convolution = 0;
    for (std::size_t k = 0; k <= j; ++k)
      convolution += first.weights()[k]*second.weights()[j-k];
    close(convolution, joined.weights()[j], 1e-11L, "LL evolution composition");
  }

  const LLNonSinglet small(.001L), twice(.002L);
  close((small.moment(3)-1)/.001L, -25.0L/12, .004L,
        "LL first-order re-expansion to P0");
  close((twice.moment(3)-2*small.moment(3)+1)/1e-6L,
        std::pow(25.0L/12,2), .02L, "LL second-order P0 tensor P0 coefficient");
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

  throws<std::invalid_argument>([&] { LLNonSinglet(-.1L, coarse); }, "backward LL time");
  throws<std::invalid_argument>([&] { evolved.moment(.5L); }, "unbounded numerical moment");
  throws<std::invalid_argument>([] { non_singlet_kernel_moment(0); }, "invalid Mellin moment");
  throws<std::invalid_argument>([] {
    LLNonSinglet(.01L,{1000,32,.1L,1e-10L});
  }, "non-power-of-two FFT grid");
  throws<std::invalid_argument>([] {
    LLNonSinglet(.01L,{1024,32,1e-6L,1e-10L});
  }, "unresolved cutoff below grid resolution");
  throws<std::invalid_argument>([] {
    LLNonSinglet(.1L,{1024,2,.01L,1e-10L});
  }, "periodic aliasing must not silently pass");
  std::cout << "LL moments, composition, expansion, endpoints and decay checks passed\n";
}
