#include "leptonff/evolution.hpp"
#include "reference_config.hpp"
#include "ll_grid_options.hpp"

#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {
using namespace leptonff;

void write_solution(const LLNonSinglet& evolution, const ReferenceInput& input) {
  const auto& d = evolution.diagnostics();
  std::cout << "{\n"
      << "      \"bins\": " << evolution.grid().bins << ",\n"
      << "      \"log_span\": " << evolution.grid().log_span << ",\n"
      << "      \"spacing\": " << evolution.grid().spacing() << ",\n"
      << "      \"soft_log_cutoff\": " << evolution.grid().soft_log_cutoff << ",\n"
      << "      \"max_wrap_probability\": " << evolution.grid().max_wrap_probability << ",\n"
      << "      \"resolved_jump_rate\": " << d.resolved_jump_rate << ",\n"
      << "      \"no_resolved_emission_weight\": " << d.no_resolved_emission_weight << ",\n"
      << "      \"wrap_probability_bound\": " << d.wrap_probability_bound << ",\n"
      << "      \"omitted_large_jump_probability_bound\": "
      << d.omitted_large_jump_probability_bound << ",\n"
      << "      \"mean_log_loss_discretisation_bound\": "
      << d.mean_log_loss_discretisation_bound << ",\n"
      << "      \"normalisation_residual\": " << d.normalisation_residual << ",\n"
      << "      \"negative_weight_mass\": " << d.negative_weight_mass << ",\n"
      << "      \"imaginary_weight_mass\": " << d.imaginary_weight_mass << ",\n"
      << "      \"acceptance\": " << decay_ll_acceptance(input.mass,input.cuts,evolution)
      << ",\n      \"moments\": [";
  for (Real N : {1.0L, 2.0L, 3.0L, 10.0L}) {
    if (N != 1) std::cout << ", ";
    std::cout << "{\"N\": " << N << ", \"numerical\": " << evolution.moment(N)
              << ", \"continuum\": " << evolution.exact_moment(N) << "}";
  }
  std::cout << "],\n      \"weights_head\": [";
  for (std::size_t j = 0; j < 16; ++j) {
    if (j != 0) std::cout << ", ";
    std::cout << evolution.weights()[j];
  }
  std::cout << "]\n    }";
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << "Usage: leptonff_ll config [bins [soft_log_cutoff [log_span]]]\n"
                   "Reports fixed-alpha non-singlet LL with grid/cutoff refinement.\n";
      return 0;
    }
    if (argc < 2 || argc > 5)
      throw std::invalid_argument("Usage: leptonff_ll config [bins [soft_log_cutoff [log_span]]]");
    ReferenceInput input;
    read_reference_config(argv[1],input);
    input.qed.validate();
    input.cuts.validate();
    const auto grid = read_ll_grid(argc,argv);
    const Real time = input.qed.expansion_parameter()
        *2*(std::log(input.qed.mu_frag)-std::log(input.qed.mu_0));
    const LLNonSinglet baseline(time,grid);
    auto fine_grid = grid;
    fine_grid.bins *= 2;
    const LLNonSinglet refined_grid(time,fine_grid);
    auto fine_cut = fine_grid;
    fine_cut.soft_log_cutoff /= 2;
    const LLNonSinglet refined_cut(time,fine_cut);
    const auto kernel = leading_kernel(Species::muon_minus,Species::muon_minus);
    const Real first_order_coefficient = 2*kernel.apply([&](Real z) {
      return isotropic_acceptance(input.mass,z,1,input.cuts);
    },0,input.integration,one_leg_cut_breakpoints(input.mass,input.cuts));
    const Real born = isotropic_acceptance(input.mass,1,1,input.cuts);

    std::cout << std::setprecision(18)
        << "{\n  \"schema_version\": 1,\n"
        << "  \"library_version\": \"" << LEPTONFF_VERSION << "\",\n"
        << "  \"build_type\": \"" << LEPTONFF_BUILD_TYPE << "\",\n"
        << "  \"compiler\": " << std::quoted(__VERSION__) << ",\n"
        << "  \"real_digits\": " << std::numeric_limits<Real>::digits << ",\n"
        << "  \"model\": \"isotropic_massless_decay_at_rest\",\n"
        << "  \"accuracy\": \"fixed_alpha_photonic_non_singlet_LL\",\n"
        << "  \"hard_Oalpha_included\": false,\n"
        << "  \"finite_Oalpha_boundary_included\": false,\n"
        << "  \"mixing_included\": false,\n"
        << "  \"units\": \"GeV_and_dimensionless_acceptances\",\n"
        << "  \"parent_mass\": " << input.mass << ",\n"
        << "  \"muon_mass\": " << input.qed.muon_mass << ",\n"
        << "  \"mu_0\": " << input.qed.mu_0 << ",\n"
        << "  \"mu_frag\": " << input.qed.mu_frag << ",\n"
        << "  \"alpha\": " << input.qed.alpha << ",\n"
        << "  \"time\": " << time << ",\n"
        << "  \"pt_min\": " << input.cuts.pt_min << ",\n"
        << "  \"eta_max\": " << input.cuts.eta_max << ",\n"
        << "  \"mass_min\": " << input.cuts.mass_min << ",\n"
        << "  \"mass_max\": " << input.cuts.mass_max << ",\n"
        << "  \"absolute_tolerance\": " << input.integration.absolute_tolerance << ",\n"
        << "  \"relative_tolerance\": " << input.integration.relative_tolerance << ",\n"
        << "  \"max_intervals\": " << input.integration.max_intervals << ",\n"
        << "  \"born_acceptance\": " << born << ",\n"
        << "  \"LL_first_order_coefficient_of_time\": " << first_order_coefficient << ",\n"
        << "  \"LL_first_order_acceptance\": " << born+time*first_order_coefficient << ",\n"
        << "  \"solutions\": {\n    \"baseline\": ";
    write_solution(baseline,input);
    std::cout << ",\n    \"refined_grid\": ";
    write_solution(refined_grid,input);
    std::cout << ",\n    \"refined_cutoff\": ";
    write_solution(refined_cut,input);
    std::cout << "\n  }\n}\n";
  } catch (const std::exception& error) {
    std::cerr << "leptonff_ll: " << error.what() << '\n';
    return 1;
  }
}
