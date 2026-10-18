#include "leptonff/fragmentation.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  using S = Species;
  QEDParameters q;
  for (Real scale : {q.muon_mass, .5L, 91.1876L}) {
    q.mu_0 = scale;
    const Real L = 2*std::log(scale/q.muon_mass);
    const auto same = initial_condition_coefficient(S::muon_minus, S::muon_minus, q);
    const auto photon = initial_condition_coefficient(S::muon_minus, S::photon, q);
    close(same.moment(1), 0, 1e-12L, "lepton number coefficient");
    close(same.moment(2), -4*L/3 - 17.0L/9, 1e-8L,
          "independent N2 lepton boundary");
    close(same.moment(3), -25*L/12 - 265.0L/72, 1e-8L,
          "independent N3 lepton boundary");
    close(photon.moment(2), 4*L/3 + 17.0L/9, 1e-8L,
          "independent N2 photon boundary");
    for (const auto parent : one_family_species) {
      Real coefficient_sum = 0, full_sum = 0, kernel_sum = 0;
      for (const auto observed : one_family_species) {
        coefficient_sum += initial_condition_coefficient(parent, observed, q).moment(2);
        full_sum += fragmentation_order1(parent, observed, q).moment(2);
        kernel_sum += leading_kernel(parent, observed).moment(2);
        const auto conjugate = [](S s) {
          return s == S::photon ? s : (s == S::muon_minus ? S::muon_plus : S::muon_minus);
        };
        close(initial_condition_coefficient(parent, observed, q).moment(3),
              initial_condition_coefficient(conjugate(parent), conjugate(observed), q).moment(3),
              1e-10L, "charge conjugation");
      }
      close(coefficient_sum, 0, 1e-8L, "parent boundary momentum");
      close(full_sum, 1, 1e-8L, "parent transported momentum");
      close(kernel_sum, 0, 1e-9L, "parent kernel momentum");
    }
    close(initial_condition_coefficient(S::photon, S::muon_plus, q).moment(2),
          L/3, 1e-10L, "photon conversion per charge");
    close(initial_condition_coefficient(S::photon, S::photon, q).moment(2),
          -2*L/3, 1e-10L, "photon survival coefficient");
    throws<std::invalid_argument>([&] { photon.moment(1); }, "soft photon number divergence");
  }

  // Starting scale cancels at O(alpha), including finite boundary pieces.
  for (const auto parent : one_family_species)
    for (const auto observed : one_family_species) {
      q.mu_0 = q.muon_mass/2;
      const Real low = fragmentation_coefficient_order1(parent, observed, q).moment(2);
      q.mu_0 = 2*q.muon_mass;
      close(fragmentation_coefficient_order1(parent, observed, q).moment(2), low,
            1e-8L, "starting scale cancellation");
      const Real original = fragmentation_coefficient_order1(parent, observed, q).moment(2);
      auto varied = q;
      varied.mu_frag *= std::exp(.25L); // delta log(mu_frag^2) = 1/2
      close(fragmentation_coefficient_order1(parent, observed, varied).moment(2)-original,
            .5L*leading_kernel(parent, observed).moment(2), 1e-8L,
            "fragmentation scale derivative");
    }
  q.alpha = 0;
  close(fragmentation_order1(S::muon_plus, S::muon_plus, q).moment(2), 1,
        1e-15L, "zero coupling identity");
  q.muon_mass = 0;
  throws<std::invalid_argument>([&] {
    initial_condition_coefficient(S::muon_plus, S::muon_plus, q);
  }, "invalid muon mass");
  std::cout << "One-family MSbar boundary and scale checks passed\n";
}
