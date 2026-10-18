#include "leptonff/hard_dressing.hpp"
#include "leptonff/vector_decay.hpp"
#include "check.hpp"
#include <numeric>

int main() {
  using namespace leptonff;
  const Real M = 91.1876L;
  const LLGrid grid{128,32,.15L,1e-10L};
  const LLNonSinglet ll(.015L,grid);
  const auto projected = evolve_hard_projection(M,M,ll);
  for (std::size_t j : {std::size_t{1},std::size_t{3},std::size_t{7},std::size_t{30}}) {
    const Real lo = std::exp(-(j+.5L)*grid.spacing());
    const Real hi = std::exp(-std::max(grid.soft_log_cutoff,(j-.5L)*grid.spacing()));
    const Real integral = vector_mass_hard_coefficient(M,M).apply(
        [=](Real z) { return z>lo && z<hi ? 1.0L : 0.0L; },0,{}, {lo,hi});
    close(projected.initial_weights[j],integral,1e-9L,"hard primitive vs distribution cell");
  }
  // Independent cyclic sum checks the FFT with a signed hard distribution.
  for (std::size_t j = 0; j < grid.bins; ++j) {
    Real direct = 0;
    for (std::size_t k = 0; k < grid.bins; ++k)
      direct += projected.initial_weights[k]*ll.weights()[(j+grid.bins-k)%grid.bins];
    close(projected.evolved_weights[j],direct,1e-10L,"signed hard evolution vs literal convolution");
  }
  close(projected.number_residual,0,1e-10L,"inclusive hard correction is 3/2");
  for (std::size_t n : {std::size_t{1},std::size_t{8},std::size_t{32}}) {
    const auto compressed = compress_ll_quadrature(ll,n);
    close(compressed.nodes[0].weight,ll.weights()[0],0,"retain LL endpoint exactly");
    close(compressed.number_residual,0,1e-12L,"compressed number moment");
    close(compressed.first_moment_residual,0,1e-12L,"compressed first z moment");
    close(compressed.second_moment_residual,0,1e-12L,"compressed second z moment");
    close(compressed.third_moment_residual,0,1e-12L,"compressed third z moment");
  }

  const FiducialCuts cuts;
  for (Real z : {.45L,.7L,.93L})
    for (Real c : {-.9L,.3L,.95L})
      close(vector_fragmented_acceptance(M,z,c,1,1,cuts),
            vector_radiative_acceptance(M,0,z,c,cuts),1e-13L,
            "fragmented measurement reduces to real event");
  for (Real z : {.55L,.8L})
    for (Real f : {.7L,.95L}) {
      close(vector_fragmented_acceptance(M,z,1,f,1,cuts),
            isotropic_acceptance(M,z*f,1,cuts),1e-13L,"first local collinear projection");
      close(vector_fragmented_acceptance(M,z,-1,f,1,cuts),
            isotropic_acceptance(M,f,z,cuts),1e-13L,"second local collinear projection");
    }

  // With only a mass cut, two LL legs collapse to U(2*t). This independent
  // one-dimensional reduction tests dressing on both muons, including delta.
  const FiducialCuts mass_only{0,100,60,120};
  const auto pair = circular_lattice_convolution(ll.weights(),ll.weights());
  const auto hard_pair = circular_lattice_convolution(projected.initial_weights,pair.weights);
  Real direct_window = 0;
  for (std::size_t j = 0; j < grid.bins; ++j) {
    const Real observed = M*std::exp(-grid.spacing()*j/2);
    if (observed>60 && observed<120) direct_window += hard_pair.weights[j];
  }
  close(decay_lattice_acceptance(M,mass_only,grid.spacing(),
         projected.evolved_weights,ll.weights()),direct_window,1e-9L,
        "mass-window hard dressing vs independent U(2t) reduction");
  const LLNonSinglet doubled(2*ll.time(),grid);
  for (std::size_t j = 0; j < grid.bins; ++j)
    close(pair.weights[j],doubled.weights()[j],1e-12L,"two-leg LL semigroup");

  QEDParameters q;
  HardDressingOptions options;
  options.ff_panels = 4;
  Real first = 0, first_dressing = 0;
  for (Real alpha : {1e-4L,2e-4L}) {
    q.alpha = alpha;
    const LLBoundaryNonSinglet u(q,grid);
    const auto fixed = fiducial_vector_order1(M,cuts,q);
    const auto match = fiducial_vector_hard_dressed_match(M,cuts,u,fixed,options);
    const Real remainder = match.matched-fixed.rate;
    const Real dressing = match.matched-match.additive_matched;
    if (alpha == 1e-4L) { first = remainder; first_dressing = dressing; }
    else {
      if (std::abs(first) < 1e-12L || std::abs(first_dressing) < 1e-12L)
        throw std::runtime_error("Unresolved hard-dressing coupling remainder");
      close(remainder/first,4,.03L,"complete dressed match starts beyond Oalpha");
      close(dressing/first_dressing,4,.03L,"hard-generated correction starts at alpha squared");
    }
  }
  q = QEDParameters{};
  q.mu_0 = q.mu_frag;
  const LLBoundaryNonSinglet no_evolution(q,grid);
  const auto fixed = fiducial_vector_order1(M,cuts,q);
  const auto match = fiducial_vector_hard_dressed_match(M,cuts,no_evolution,fixed,options);
  close(match.projected_dressing_shift,0,1e-11L,"identity hard projection evolution");
  close(match.geometry_dressing_shift,0,1e-13L,"identity hard geometry evolution");
  close(match.matched,fixed.rate,1e-12L,"dressed no-evolution match equals full Oalpha");

  q = QEDParameters{};
  const LLBoundaryNonSinglet physical(q,grid);
  const auto inclusive_fixed = fiducial_vector_order1(M,{0,100,0,120},q);
  const auto inclusive = fiducial_vector_hard_dressed_match(
      M,{0,100,0,120},physical,inclusive_fixed,options);
  close(inclusive.dressed_hard_coefficient,1.5L,1e-9L,"inclusive dressed virtual and plus normalisation");
  close(inclusive.matched,1+q.expansion_parameter()*1.5L,1e-10L,
        "inclusive LL/boundary number conservation with hard matching");
  throws<std::invalid_argument>([&] {
    fiducial_vector_hard_dressed_match(M,cuts,physical,fixed,options);
  },"reject incompatible cached hard coefficient");
  throws<std::invalid_argument>([&] { compress_ll_quadrature(ll,0); },"invalid FF quadrature");
  throws<std::invalid_argument>([&] {
    circular_lattice_convolution({1,0},{1,0});
  },"invalid cyclic convolution size");
  std::cout << "Hard LL dressing, signed convolution, moments and expansion checks passed\n";
}
