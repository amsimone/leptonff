#include "leptonff/nlo_qed.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  using S = Species;
  for (S parent : one_family_species) {
    Real momentum = 0;
    for (S observed : one_family_species) {
      const auto kernel = next_to_leading_kernel(parent,observed);
      momentum += kernel.moment(2);
      auto conjugate = [](S s) {
        return s == S::photon ? s : s == S::muon_minus ? S::muon_plus : S::muon_minus;
      };
      close(kernel.moment(3),
            next_to_leading_kernel(conjugate(parent),conjugate(observed)).moment(3),
            1e-12L,"full matrix charge conjugation");
    }
    close(momentum,0,2e-8L,"NLO per-parent momentum conservation");
  }
  const auto photon = next_to_leading_kernel(S::photon,S::photon);
  close(photon.moment(2),-65.0L/27,2e-9L,"analytic photon N2 including delta");
  close(next_to_leading_kernel(S::photon,S::muon_minus).moment(2),65.0L/54,
        2e-9L,"analytic conversion N2 per charge");
  // Analytic regular photon moments, evaluated without any dilogarithm.
  for (Real N : {2.L,3.L,5.L}) {
    const Real expected = -1-4/N+12/(N+1)-164/(9*(N+2))+92/(9*(N-1))
        -10/(N*N)-14/std::pow(N+1,2)-16/(3*std::pow(N+2,2))
        -16/(3*std::pow(N-1,2))+4/std::pow(N,3)+4/std::pow(N+1,3);
    close(photon.moment(N),expected,2e-9L,"analytic photon Mellin moment");
  }
  throws<std::invalid_argument>([&] { photon.moment(1); },"photon NLO number divergence");
  throws<std::invalid_argument>([] {
    next_to_leading_kernel(static_cast<S>(99),S::photon);
  },"unknown NLO species");
  std::cout << "Complete one-family time-like NLO matrix checks passed\n";
}
