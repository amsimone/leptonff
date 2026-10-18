#include "leptonff/vector_decay.hpp"
#include "reference_config.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace leptonff;

namespace {
void print_result(const VectorMassResult& r) {
  std::cout << "{\"born\":" << r.born
      << ",\"hard_coefficient\":" << r.hard_coefficient
      << ",\"ff_coefficient\":" << r.ff_coefficient
      << ",\"matched_coefficient\":" << r.matched_coefficient
      << ",\"leading_power_coefficient\":" << r.leading_power_coefficient
      << ",\"leading_power_rate\":" << r.leading_power_rate
      << ",\"massive_born\":" << r.massive_born
      << ",\"massive_coefficient\":" << r.massive_coefficient
      << ",\"massive_rate\":" << r.massive_rate << '}';
}
} // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << "Usage: leptonff_hard [vector_mass.conf]\n"
                   "Angular-inclusive vector-current mass-window benchmark; no pT/eta cuts.\n";
      return 0;
    }
    if (argc > 2) throw std::invalid_argument("Usage: leptonff_hard [vector_mass.conf]");
    ReferenceInput input;
    input.cuts.pt_min = 0;
    input.cuts.eta_max = 100;
    if (argc == 2) read_reference_config(argv[1],input);
    input.cuts.validate();
    input.qed.validate();
    if (input.cuts.pt_min != 0 || std::tanh(input.cuts.eta_max) != 1)
      throw std::invalid_argument("Mass-only benchmark requires pt_min=0 and an inactive eta cut (use eta_max=100)");
    const MassWindow window{input.cuts.mass_min,input.cuts.mass_max};
    const auto result = vector_mass_window(input.mass,window,input.qed,input.integration);
    auto tighter = input.integration;
    tighter.absolute_tolerance /= 10;
    tighter.relative_tolerance /= 10;
    const auto refined = vector_mass_window(input.mass,window,input.qed,tighter);
    std::cout << std::setprecision(18)
        << "{\n  \"schema_version\":1,\n  \"library_version\":\"" << LEPTONFF_VERSION
        << "\",\n  \"build_type\":\"" << LEPTONFF_BUILD_TYPE
        << "\",\n  \"model\":\"angular_inclusive_unpolarised_vector_current\",\n"
        << "  \"accuracy\":\"Oalpha_hard_plus_FF_leading_power_with_massive_reference\",\n"
        << "  \"normalisation\":\"massless_vector_Born\",\n"
        << "  \"coefficient_expansion\":\"alpha_over_2pi\",\n"
        << "  \"scheme\":\"MSbar\",\n  \"coupling\":\"fixed_alpha\",\n"
        << "  \"hard_Oalpha_included\":true,\n  \"resummed\":false,\n"
        << "  \"angular_cuts\":false,\n  \"axial_current\":false,\n"
        << "  \"units\":\"GeV_and_rates_relative_to_massless_vector_Born\",\n"
        << "  \"compiler\":" << std::quoted(__VERSION__)
        << ",\n  \"real_digits\":" << std::numeric_limits<Real>::digits
        << ",\n  \"parent_mass\":" << input.mass
        << ",\n  \"muon_mass\":" << input.qed.muon_mass
        << ",\n  \"mu_0\":" << input.qed.mu_0
        << ",\n  \"mu_frag\":" << input.qed.mu_frag
        << ",\n  \"alpha\":" << input.qed.alpha
        << ",\n  \"pt_min\":" << input.cuts.pt_min
        << ",\n  \"eta_max\":" << input.cuts.eta_max
        << ",\n  \"mass_min\":" << window.lower
        << ",\n  \"mass_max\":" << window.upper
        << ",\n  \"absolute_tolerance\":" << input.integration.absolute_tolerance
        << ",\n  \"relative_tolerance\":" << input.integration.relative_tolerance
        << ",\n  \"max_intervals\":" << input.integration.max_intervals
        << ",\n  \"result\":";
    print_result(result);
    std::cout << ",\n  \"leading_power_refinement_shift\":"
        << refined.leading_power_rate-result.leading_power_rate
        << ",\n  \"massive_refinement_shift\":" << refined.massive_rate-result.massive_rate
        << ",\n  \"scale_scan\":[";
    bool first = true;
    for (Real factor : {.5L,1.0L,2.0L}) {
      auto q = input.qed;
      q.mu_frag = factor*input.mass;
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"mu_frag\":" << q.mu_frag << ",\"result\":";
      print_result(vector_mass_window(input.mass,window,q,input.integration));
      std::cout << '}';
    }
    std::cout << "],\n  \"mass_scan\":[";
    first = true;
    for (Real fraction : {1.0L,.5L,.25L,.125L}) {
      auto q = input.qed;
      q.muon_mass *= fraction;
      q.mu_0 *= fraction;
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"muon_mass\":" << q.muon_mass << ",\"result\":";
      print_result(vector_mass_window(input.mass,window,q,input.integration));
      std::cout << '}';
    }
    std::cout << "]\n}\n";
  } catch (const std::exception& error) {
    std::cerr << "leptonff_hard: " << error.what() << '\n';
    return 1;
  }
}
