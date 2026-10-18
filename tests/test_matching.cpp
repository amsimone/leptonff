#include "leptonff/decay.hpp"
#include "leptonff/evolution.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  using S = Species;
  QEDParameters q;
  const LLGrid coarse{512,32,.04L,1e-10L};
  const LLBoundaryNonSinglet fixture(q,coarse);
  const auto& u = fixture.evolution();
  const auto& b = fixture.initial_boundary_weights();
  const auto& evolved = fixture.evolved_boundary_weights();
  for (std::size_t j : {std::size_t{1},std::size_t{7},std::size_t{18}}) {
    const Real lo = std::exp(-(j+.5L)*coarse.spacing());
    const Real hi = std::exp(-std::max(coarse.soft_log_cutoff,(j-.5L)*coarse.spacing()));
    const Real independent = initial_condition_coefficient(S::muon_minus,S::muon_minus,q)
        .apply([=](Real z) { return z>lo && z<hi ? 1.0L : 0.0L; },0,{}, {lo,hi});
    close(b[j],independent,1e-9L,"analytic boundary cell vs plus-distribution action");
  }
  for (std::size_t j = 0; j < 64; ++j) {
    Real direct = 0;
    for (std::size_t k = 0; k <= j; ++k) direct += u.weights()[k]*b[j-k];
    close(evolved[j],direct,1e-10L,"signed boundary evolution vs direct convolution");
  }
  close(fixture.boundary_moment(1),0,1e-10L,"boundary lepton number");
  close(fixture.apply([](Real) { return 1; }),1,1e-12L,"LL plus boundary lepton number");

  const LLBoundaryNonSinglet fine(q);
  const LLBoundaryNonSinglet refined(q,{524288,32,1e-4L,1e-10L});
  for (Real N : {2.0L,3.0L}) {
    const Real exact = fine.exact_boundary_moment(N);
    if (!(std::abs(refined.boundary_moment(N)-exact) <
          std::abs(fine.boundary_moment(N)-exact)))
      throw std::runtime_error("Boundary moment failed regulator convergence");
    close(refined.boundary_moment(N),exact,.025L,"boundary evolved moment reference");
  }

  const FiducialCuts cuts;
  Real direct = 0;
  for (std::size_t i = 0; i < coarse.bins; ++i)
    for (std::size_t j = 0; j < coarse.bins; ++j)
      direct += u.weights()[i]*evolved[j]
          *isotropic_acceptance(91.1876L,u.z(i),u.z(j),cuts);
  close(decay_lattice_acceptance(91.1876L,cuts,coarse.spacing(),u.weights(),evolved),
        direct,1e-10L,"signed bilinear two-leg measurement");

  // Even at coarse resolution the additive match removes its complete O(a)
  // overlap. The remaining difference from fixed order must start at a^2.
  q.alpha = 1e-4L;
  const auto small = decay_ff_ll_match(91.1876L,cuts,LLBoundaryNonSinglet(q,coarse));
  Real remainder1 = small.matched-small.born
      -q.expansion_parameter()*small.fixed_order_coefficient;
  q.alpha *= 2;
  const auto twice = decay_ff_ll_match(91.1876L,cuts,LLBoundaryNonSinglet(q,coarse));
  Real remainder2 = twice.matched-twice.born
      -q.expansion_parameter()*twice.fixed_order_coefficient;
  if (std::abs(remainder1) < 1e-12L)
    throw std::runtime_error("Matching expansion test has no resolved remainder");
  close(remainder2/remainder1,4,.02L,"matched FF fixed-order overlap starts at a^2");

  q = QEDParameters{};
  q.mu_0 = q.mu_frag; // no evolution: only a single boundary insertion remains
  const LLBoundaryNonSinglet equal_scales(q,coarse);
  const auto equal_match = decay_ff_ll_match(91.1876L,cuts,equal_scales);
  close(equal_match.matched,equal_match.born
        +q.expansion_parameter()*equal_match.fixed_order_coefficient,1e-12L,
        "no-evolution match equals full Oalpha FF");
  auto full_weights = equal_scales.evolution().weights();
  for (std::size_t j = 0; j < coarse.bins; ++j)
    full_weights[j] += q.expansion_parameter()*equal_scales.evolved_boundary_weights()[j];
  const Real uncontrolled_product = decay_lattice_acceptance(
      91.1876L,cuts,coarse.spacing(),full_weights,full_weights);
  const Real unwanted = std::pow(q.expansion_parameter(),2)*decay_lattice_acceptance(
      91.1876L,cuts,coarse.spacing(),equal_scales.evolved_boundary_weights(),
      equal_scales.evolved_boundary_weights());
  close(uncontrolled_product-equal_match.resummed_with_boundary,unwanted,1e-12L,
        "global single-boundary truncation excludes a^2 d1 x d1");
  if (std::abs(unwanted) < 1e-8L)
    throw std::runtime_error("Boundary double-insertion check is not resolved");

  q.alpha = 0;
  close(decay_ff_ll_match(91.1876L,cuts,LLBoundaryNonSinglet(q,coarse)).matched,
        isotropic_acceptance(91.1876L,1,1,cuts),1e-12L,"matched zero-coupling limit");
  throws<std::invalid_argument>([&] {
    decay_lattice_acceptance(91.1876L,cuts,.1L,{1,0},{1});
  },"incompatible signed lattices");
  std::cout << "Evolved signed boundary and FF-only additive matching checks passed\n";
}
