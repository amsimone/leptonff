#include "leptonff/pair_moment.hpp"
#include "leptonff/vector_decay.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const Real M = 91.1876L;
  MomentEvolutionParameters p;
  MomentEvolutionOptions options;
  options.step_tolerance = 1e-13L;
  options.kernel_integration = {1e-10L,1e-10L,8192};
  for (auto mode : {QEDCouplingMode::fixed,QEDCouplingMode::one_loop,QEDCouplingMode::two_loop}) {
    p.coupling = mode;
    for (Real N : {2.0L,3.0L,5.0L}) {
      const auto central = vector_pair_moment(M,N,p,options);
      close(central.coefficient1,vector_mass_leading_power_coefficient(M,p.muon_mass)
          .moment(N,options.kernel_integration),1e-8L,"complete Oalpha pair-mass coefficient");
      auto q = p;
      q.mu_frag = 2*M;
      const auto varied = vector_pair_moment(M,N,q,options);
      close(varied.coefficient1,central.coefficient1,1e-8L,"pair-moment Oalpha scale cancellation");
      close(varied.a_hard,central.a_hard,1e-15L,"hard coupling fixed at M across mu_frag scan");
      // At alpha^2 the residual mu_frag dependence may have constants, but
      // cannot contain log(M/m). All NLL mass logs must cancel in the shift.
      const Real shift = varied.coefficient2-central.coefficient2;
      q.muon_mass /= 2; q.mu_0 = q.muon_mass;
      const auto low_varied = vector_pair_moment(M,N,q,options);
      q.mu_frag = M;
      const auto low_central = vector_pair_moment(M,N,q,options);
      close(low_varied.coefficient2-low_central.coefficient2,shift,2e-7L,
          "factorisation variation has no retained alpha^2 mass logarithm");
      p.alpha_msbar_at_mu0 = .001L;
      const auto small = vector_pair_moment(M,N,p,options);
      const Real remainder = small.matched-1-small.a_initial*small.coefficient1
          -small.a_initial*small.a_initial*small.coefficient2;
      p.alpha_msbar_at_mu0 *= 2;
      const auto doubled = vector_pair_moment(M,N,p,options);
      const Real remainder2 = doubled.matched-1-doubled.a_initial*doubled.coefficient1
          -doubled.a_initial*doubled.a_initial*doubled.coefficient2;
      close(remainder2/remainder,8,.2L,"matched pair moment cubic expansion remainder");
      p.alpha_msbar_at_mu0 = 1/137.035999084L;
    }
  }
  p.alpha_msbar_at_mu0 = 0;
  close(vector_pair_moment(M,2,p,options).matched,1,0,"zero-coupling pair moment");
  throws<std::invalid_argument>([&] { vector_pair_moment(M,1,p); },"pair count unresolved at N1");
  std::cout << "Matched pair moments, hard coupling, NLL channel and expansion checks passed\n";
}
