#include "leptonff/nlo_lepton.hpp"
#include "check.hpp"

namespace {
leptonff::Real harmonic(unsigned n, unsigned power) {
  leptonff::Real result = 0;
  for (unsigned j = 1; j <= n; ++j) result += 1/std::pow(leptonff::Real(j), power);
  return result;
}
}

int main() {
  using namespace leptonff;
  using C = NLOLeptonComponent;
  using S = Species;
  const Real zeta2 = std::pow(std::acos(-1.0L),2)/6;
  const Real k = 13.0L/4-3*zeta2+2*1.2020569031595942854L;
  // Source eqs. (5.16)--(5.21): conservation applies to the combination,
  // NOT to the source's photonic diagram class alone.
  close(nlo_lepton_component(C::photonic).moment(1), -k, 2e-9L, "photonic N1");
  close(nlo_lepton_component(C::nonsinglet_pair).moment(1), 0, 2e-9L, "NS pair N1");
  close(nlo_lepton_component(C::same_charge_interference).moment(1), 2*k,
        2e-9L, "same interference N1");
  close(nlo_lepton_component(C::opposite_charge_interference).moment(1), k,
        2e-9L, "opposite interference N1");
  close(nlo_valence_kernel().moment(1), 0, 2e-9L, "full valence charge conservation");

  for (unsigned n : {2u, 3u, 5u}) {
    const Real N = n;
    const Real ns = 20*harmonic(n-1,1)/9
        -2*(harmonic(n-1,2)+harmonic(n+1,2))/3-1.0L/6
        -2/(9*N)+22/(9*(N+1));
    const Real ps = 2/(N*N*N)+2/std::pow(N+1,3)+5/(N*N)
        +9/std::pow(N+1,2)+8/(3*std::pow(N+2,2))
        -8/N-20/(9*(N-1))+4/(N+1)+56/(9*(N+2));
    close(nlo_lepton_component(C::nonsinglet_pair).moment(N), ns, 2e-9L,
          "analytic NS pair moment");
    close(nlo_lepton_component(C::singlet_pair).moment(N), ps, 2e-9L,
          "analytic singlet pair moment");
    const auto same = nlo_lepton_kernel(S::muon_minus, S::muon_minus);
    const auto opposite = nlo_lepton_kernel(S::muon_minus, S::muon_plus);
    close(same.moment(N)-opposite.moment(N), nlo_valence_kernel().moment(N),
          3e-9L, "singlet cancellation in valence");
    close(same.moment(N), nlo_lepton_kernel(S::muon_plus,S::muon_plus).moment(N),
          1e-12L, "same-charge conjugation");
    close(opposite.moment(N), nlo_lepton_kernel(S::muon_plus,S::muon_minus).moment(N),
          1e-12L, "opposite-charge conjugation");
  }

  const auto decomposed = nlo_lepton_component(C::photonic)
      +nlo_lepton_component(C::nonsinglet_pair)
      +nlo_lepton_component(C::singlet_pair)
      +nlo_lepton_component(C::same_charge_interference);
  for (Real lower : {0.0L, .37L, .91L})
    close(nlo_lepton_kernel(S::muon_minus,S::muon_minus).apply(
              [](Real z) { return z*z; }, lower),
          decomposed.apply([](Real z) { return z*z; }, lower), 2e-9L,
          "independent log-only assembly at finite bounds");
  throws<std::invalid_argument>([] {
    nlo_lepton_component(C::singlet_pair).moment(1);
  }, "singlet number divergence");
  throws<std::invalid_argument>([] {
    nlo_lepton_kernel(S::muon_minus,S::muon_minus).moment(1);
  }, "inclusive lepton number divergence");
  throws<std::invalid_argument>([] {
    nlo_lepton_kernel(S::muon_minus,S::muon_plus).moment(1);
  }, "inclusive antilepton number divergence");
  for (S parent : one_family_species)
    for (S observed : one_family_species)
      if (parent == S::photon || observed == S::photon)
        throws<std::invalid_argument>([=] { nlo_lepton_kernel(parent,observed); },
                                      "unimplemented photon kernel");
  throws<std::invalid_argument>([] {
    nlo_lepton_component(static_cast<C>(99));
  }, "unknown NLO component");
  std::cout << "Time-like NLO lepton block and valence checks passed\n";
}
