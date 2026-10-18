#include "leptonff/decay.hpp"
#include "check.hpp"
#include <algorithm>
#include <limits>

int main() {
  using namespace leptonff;
  constexpr Real M = 91.1876L;
  const FiducialCuts cuts;
  for (Real c : {-.9L, -.3L, 0.0L, .2L, .95L}) {
    const auto pair = born_decay(M, c, .8L);
    close(invariant_mass(pair.first, pair.second), M, 1e-12L, "Born pair mass");
    for (Real z1 : {.25L, .7L, 1.0L}) {
      const Real z2 = .83L;
      const auto first = pair.first.scaled(z1);
      const auto second = pair.second.scaled(z2);
      close(first.pt(), z1*pair.first.pt(), 1e-12L, "pt rescaling");
      close(first.eta(), pair.first.eta(), 1e-12L, "eta invariance");
      close(invariant_mass(first, second), M*std::sqrt(z1*z2), 1e-12L,
            "two-leg mass rescaling");
    }
  }

  const Real angular_bound = std::sqrt(1-std::pow(2*cuts.pt_min/M,2));
  const Real angular = integrate([&](Real c) {
    const auto pair = born_decay(M, c);
    return cuts.accepts(pair.first, pair.second) ? .5L : 0.0L;
  }, -1, 1, {}, {-angular_bound, angular_bound}).value;
  close(isotropic_acceptance(M,1,1,cuts), angular, 1e-12L,
        "analytic acceptance vs independent angular cut evaluation");
  const Distribution identity(1);
  close(apply_pair(identity, identity, [&](Real z1, Real z2) {
    return isotropic_acceptance(M,z1,z2,cuts);
  }), angular, 1e-12L, "identity pair convolution");

  auto eta_cuts = cuts;
  eta_cuts.pt_min = 0;
  eta_cuts.eta_max = .4L;
  close(isotropic_acceptance(M,1,1,eta_cuts), std::tanh(.4L), 1e-15L,
        "eta-limited acceptance");
  close(isotropic_acceptance(M,.3L,1,cuts), 0, 1e-15L, "migration out of cuts");
  // A parent above the mass window can migrate into it after fragmentation.
  close(isotropic_acceptance(130,1,1,cuts), 0, 1e-15L, "parent outside mass window");
  if (!(isotropic_acceptance(130,.7L,1,cuts) > 0))
    throw std::runtime_error("lost migration into fiducial window");

  QEDParameters qed;
  const Real coefficient = decay_ff_coefficient_order1(M, cuts, qed);
  if (!(coefficient < 0)) throw std::runtime_error("expected radiative loss in toy cuts");
  const auto d1 = fragmentation_coefficient_order1(
      Species::muon_minus,Species::muon_minus,qed);
  const auto points = one_leg_cut_breakpoints(M,cuts);
  close(apply_pair_order1(d1,d1,qed.expansion_parameter(), [&](Real z1, Real z2) {
    return isotropic_acceptance(M,z1,z2,cuts);
  }, {}, points, points), angular+qed.expansion_parameter()*coefficient, 1e-10L,
        "explicitly expanded two-leg decay");
  auto inclusive = cuts;
  inclusive.pt_min = 0;
  inclusive.mass_min = 0;
  inclusive.eta_max = 20;
  close(decay_ff_coefficient_order1(M,inclusive,qed), 0, 1e-12L,
        "inclusive lepton number cancellation");
  auto tighter = IntegrationOptions{};
  tighter.absolute_tolerance /= 10;
  tighter.relative_tolerance /= 10;
  close(decay_ff_coefficient_order1(M,cuts,qed,tighter), coefficient, 1e-7L,
        "fiducial coefficient numerical refinement");
  // Literal signed two-leg cell integration, independent of the efficient
  // prefix/transition algorithm. Includes an infinite upper loss bound.
  const std::vector<Real> a{.7L,-.1L,.05L,.2L,-.04L,.1L};
  const std::vector<Real> b{1,.03L,-.06L,.01L,.03L,-.02L};
  const Real h=.13L;
  for(FiducialCuts window:{FiducialCuts{0,100,0,120},FiducialCuts{20,1.5L,68,88},
                           FiducialCuts{0,.7L,75,120}}) {
    Real literal=0;
    auto angular_cuts=window; angular_cuts.mass_min=0; angular_cuts.mass_max=2*M;
    const Real low=2*std::log(M/window.mass_max);
    const Real high=window.mass_min==0 ? std::numeric_limits<Real>::infinity()
        : 2*std::log(M/window.mass_min);
    for(std::size_t i=0;i<a.size();++i) for(std::size_t j=0;j<b.size();++j) {
      const Real f=std::exp(-h*i),g=std::exp(-h*j),u=h*(i+j);
      const Real fraction=i+j==0 ? (window.mass_min<M && M<window.mass_max ? 1 : 0)
          : std::max(Real{0},std::min(u+h/2,high)-std::max(u-h/2,low))/h;
      literal+=a[i]*b[j]*fraction*isotropic_acceptance(M,f,g,angular_cuts);
    }
    close(decay_lattice_acceptance(M,window,h,a,b,true),literal,1e-13L,
          "conservative mass-cell integration vs literal signed tensor");
  }
  throws<std::invalid_argument>([] { born_decay(M, 1.1L); }, "invalid cosine");
  throws<std::invalid_argument>([&] { isotropic_acceptance(M,-.1L,1,cuts); },
                                "invalid fragmentation fraction");
  std::cout << "Born kinematics and decay-convolution checks passed\n";
}
