#include "leptonff/fiducial_vector.hpp"
#include "leptonff/evolution.hpp"
#include "leptonff/vector_current.hpp"
#include "leptonff/vector_decay.hpp"
#include "check.hpp"
#include <algorithm>

int main() {
  using namespace leptonff;
  const Real M = 91.1876L, pi = std::acos(-1.0L);
  for (Real b1 : {.1L,.4L,.8L,.99L})
    for (Real b2 : {.2L,.5L,.9L})
      for (Real angle : {.0001L,.1L,.7L,1.4L}) {
        const Real cosine = std::cos(angle), sine = std::sin(angle);
        const Real independent = integrate([&](Real t) {
          const Real v = std::sqrt(1-t*t)*sine;
          const Real lo = std::clamp((-b2-t*cosine)/v,Real{-1},Real{1});
          const Real hi = std::clamp((b2-t*cosine)/v,Real{-1},Real{1});
          return (std::acos(lo)-std::acos(hi))/pi;
        },0,b1,{1e-8L,1e-8L,8192},
        {std::abs(b2*cosine-std::sqrt(1-b2*b2)*sine),
         b2*cosine+std::sqrt(1-b2*b2)*sine}).value;
        close(rotational_band_acceptance(b1,b2,cosine),independent,1e-7L,
              "spherical caps vs independent azimuth integration");
      }
  for (Real m : {0.0L,.1056583755L,2.0L,25.0L})
    for (Real z : {.45L,.7L,.93L})
      for (Real c : {-.9L,.3L,.95L})
        close(vector_scalar_real_density(M,m,z,c),vector_real_density(M,m,z,c),
              1e-8L,"scalar current vs independent Dirac trace, massive and massless");
  const FiducialCuts cuts;
  for (Real z : {.35L,.55L,.7L,.9L})
    close(vector_radiative_acceptance(M,0,z,1,cuts),
          isotropic_acceptance(M,z,1,cuts),1e-14L,"local collinear projection");
  QEDParameters q;
  HardDecayOptions options;
  const auto fixed = fiducial_vector_order1(M,cuts,q,options);
  options.soft_cutoff /= 2;
  options.collinear_cutoff /= 2;
  const auto refined = fiducial_vector_order1(M,cuts,q,options);
  close(fixed.rate,refined.rate,2e-8L,"finite remainder regulator refinement");
  options.integration.absolute_tolerance /= 4;
  options.integration.relative_tolerance /= 4;
  close(fiducial_vector_order1(M,cuts,q,options).rate,refined.rate,2e-8L,
        "nested quadrature refinement");
  const auto central_tight = fiducial_vector_order1(M,cuts,q,options);
  for (Real factor : {.5L,2.0L}) {
    q.mu_frag = factor*M;
    const auto shifted = fiducial_vector_order1(M,cuts,q,options);
    close(shifted.full_coefficient,central_tight.full_coefficient,1e-7L,
          "fiducial hard plus FF scale cancellation");
  }
  q = QEDParameters{};
  for (MassWindow window : {MassWindow{60,120},MassWindow{70,85},MassWindow{0,120}}) {
    const FiducialCuts inclusive{0,100,window.lower,window.upper};
    const auto mass_only = vector_mass_window(M,window,q);
    const auto projected = fiducial_vector_order1(M,inclusive,q);
    close(projected.rate,mass_only.leading_power_rate,1e-10L,
          "angular-inclusive leading-power reduction");
    close(massive_fiducial_vector_order1(M,inclusive,q).rate,mass_only.massive_rate,
          1e-9L,"angular-inclusive massive reduction");
  }
  // Tied eta thresholds produce linear mass sensitivity and sqrt(delta)
  // collinear-sliver errors. Test with a much smaller angular sliver rather
  // than imposing the pT fixture's quadratic mass envelope.
  auto alternate_options = options;
  alternate_options.soft_cutoff = 1e-7L;
  alternate_options.collinear_cutoff = 1e-10L;
  for (FiducialCuts alternate : {FiducialCuts{0,.8L,60,120},
                                 FiducialCuts{25,2.4L,70,85}}) {
    Real previous_alternate = 1;
    for (Real factor : {1.0L,.5L,.25L}) {
      auto alternate_q = q;
      alternate_q.muon_mass *= factor;
      alternate_q.mu_0 *= factor;
      const auto lp = fiducial_vector_order1(M,alternate,alternate_q,alternate_options);
      const auto massive = massive_fiducial_vector_order1(M,alternate,alternate_q,alternate_options);
      const Real difference = std::abs(lp.rate-massive.rate);
      if (difference > .65L*previous_alternate || difference > 1e-5L)
        throw std::runtime_error("Eta-only or real-tail massive limit fails convergence");
      previous_alternate = difference;
      if (alternate.pt_min == 0 && lp.geometry_coefficient >= 0)
        throw std::runtime_error("Separated directions must reduce eta-only acceptance");
    }
  }
  // A genuine mass scan at the same cuts; use a smaller soft sliver so its
  // truncation does not mask the power-suppressed differences.
  options.soft_cutoff = 1e-7L;
  Real previous = 1;
  for (Real factor : {1.0L,.5L,.25L}) {
    q.muon_mass = .1056583755L*factor;
    q.mu_0 = q.muon_mass;
    const auto lp = fiducial_vector_order1(M,cuts,q,options);
    const auto massive = massive_fiducial_vector_order1(M,cuts,q,options);
    const Real difference = std::abs(massive.rate-lp.rate);
    if (difference > .6L*previous || difference > 5e-6L)
      throw std::runtime_error("Fiducial massive reference fails leading-power convergence");
    previous = difference;
  }
  // The combined match must remove the complete lattice O(alpha) overlap.
  const LLGrid grid{512,32,.04L,1e-10L};
  q = QEDParameters{};
  Real first_remainder = 0;
  for (Real alpha : {1e-4L,2e-4L}) {
    q.alpha = alpha;
    const auto fo = fiducial_vector_order1(M,cuts,q);
    const auto match = fiducial_vector_ll_match(M,cuts,LLBoundaryNonSinglet(q,grid),fo);
    const Real remainder = match.matched-fo.rate;
    if (alpha == 1e-4L) first_remainder = remainder;
    else close(remainder/first_remainder,4,.02L,"full hard+FF overlap starts at alpha squared");
  }
  q = QEDParameters{};
  q.mu_0 = q.mu_frag;
  const auto no_evolution = fiducial_vector_order1(M,cuts,q);
  close(fiducial_vector_ll_match(M,cuts,LLBoundaryNonSinglet(q,grid),no_evolution).matched,
        no_evolution.rate,1e-12L,"no-evolution match equals full fixed order");
  throws<std::invalid_argument>([&] {
    fiducial_vector_ll_match(M,cuts,LLBoundaryNonSinglet(q,grid),fixed);
  },"reject stale hard coefficient configuration");
  q.alpha = 0;
  const auto zero = fiducial_vector_order1(M,cuts,q);
  close(fiducial_vector_ll_match(M,cuts,LLBoundaryNonSinglet(q,grid),zero).matched,
        zero.born,1e-12L,"zero-coupling identity");
  throws<std::invalid_argument>([&] {
    fiducial_vector_order1(M,{25,2.4L,60,M},q);
  },"unsafe upper Born endpoint cut");
  throws<std::invalid_argument>([&] { vector_real_angular_integral(M,0,.7L); },
                                "massless angular divergence");
  std::cout << "Fiducial vector hard matching, massive and NLO+LL checks passed\n";
}
