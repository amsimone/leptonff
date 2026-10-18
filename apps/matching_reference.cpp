#include "leptonff/evolution.hpp"
#include "reference_config.hpp"
#include "ll_grid_options.hpp"

#include <iomanip>
#include <iostream>
#include <limits>

namespace {
using namespace leptonff;

void write_head(const std::vector<Real>& weights) {
  std::cout << '[';
  for (std::size_t j = 0; j < 16; ++j) {
    if (j != 0) std::cout << ", ";
    std::cout << weights[j];
  }
  std::cout << ']';
}

void write_solution(const LLBoundaryNonSinglet& boundary, const ReferenceInput& input) {
  const auto& ll = boundary.evolution();
  const auto& d = boundary.diagnostics();
  const auto match = decay_ff_ll_match(input.mass,input.cuts,boundary,input.integration);
  std::cout << "{\n"
      << "      \"bins\": " << ll.grid().bins << ",\n"
      << "      \"log_span\": " << ll.grid().log_span << ",\n"
      << "      \"spacing\": " << ll.grid().spacing() << ",\n"
      << "      \"soft_log_cutoff\": " << ll.grid().soft_log_cutoff << ",\n"
      << "      \"max_wrap_probability\": " << ll.grid().max_wrap_probability << ",\n"
      << "      \"resolved_jump_rate\": " << ll.diagnostics().resolved_jump_rate << ",\n"
      << "      \"LL_wrap_probability_bound\": " << ll.diagnostics().wrap_probability_bound << ",\n"
      << "      \"boundary_wrap_coefficient_bound\": " << d.wrap_coefficient_bound << ",\n"
      << "      \"boundary_omitted_tail_coefficient_bound\": " << d.omitted_tail_coefficient_bound << ",\n"
      << "      \"boundary_initial_number_residual\": " << d.initial_number_residual << ",\n"
      << "      \"boundary_evolved_number_residual\": " << d.evolved_number_residual << ",\n"
      << "      \"boundary_imaginary_weight_mass\": " << d.imaginary_weight_mass << ",\n"
      << "      \"acceptance\": " << match.ll << ",\n"
      << "      \"evolved_boundary_coefficient\": " << match.evolved_boundary_coefficient << ",\n"
      << "      \"resummed_with_boundary\": " << match.resummed_with_boundary << ",\n"
      << "      \"lattice_overlap_coefficient\": " << match.lattice_overlap_coefficient << ",\n"
      << "      \"fixed_order_coefficient\": " << match.fixed_order_coefficient << ",\n"
      << "      \"additive_correction\": " << match.additive_correction << ",\n"
      << "      \"matched_FF_acceptance\": " << match.matched << ",\n"
      << "      \"boundary_moments\": [";
  for (Real N : {1.0L,2.0L,3.0L}) {
    if (N != 1) std::cout << ", ";
    std::cout << "{\"N\": " << N << ", \"numerical\": " << boundary.boundary_moment(N)
              << ", \"continuum\": " << boundary.exact_boundary_moment(N,input.integration) << '}';
  }
  std::cout << "],\n      \"weights_head\": ";
  write_head(ll.weights());
  std::cout << ",\n      \"initial_boundary_head\": ";
  write_head(boundary.initial_boundary_weights());
  std::cout << ",\n      \"evolved_boundary_head\": ";
  write_head(boundary.evolved_boundary_weights());
  std::cout << "\n    }";
}
} // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << "Usage: leptonff_match config [bins [soft_log_cutoff [log_span]]]\n"
                   "LL plus one evolved Oalpha boundary, matched to Oalpha FF only.\n";
      return 0;
    }
    if (argc < 2 || argc > 5)
      throw std::invalid_argument("Usage: leptonff_match config [bins [soft_log_cutoff [log_span]]]");
    ReferenceInput input;
    read_reference_config(argv[1],input);
    input.qed.validate();
    input.cuts.validate();
    const auto grid = read_ll_grid(argc,argv);
    const LLBoundaryNonSinglet baseline(input.qed,grid);
    auto fine_grid = grid;
    fine_grid.bins *= 2;
    const LLBoundaryNonSinglet refined_grid(input.qed,fine_grid);
    auto fine_cut = fine_grid;
    fine_cut.soft_log_cutoff /= 2;
    const LLBoundaryNonSinglet refined_cut(input.qed,fine_cut);
    const Real born = isotropic_acceptance(input.mass,1,1,input.cuts);
    const Real fixed = decay_ff_coefficient_order1(input.mass,input.cuts,input.qed,input.integration);

    std::cout << std::setprecision(18)
        << "{\n  \"schema_version\": 1,\n"
        << "  \"library_version\": \"" << LEPTONFF_VERSION << "\",\n"
        << "  \"build_type\": \"" << LEPTONFF_BUILD_TYPE << "\",\n"
        << "  \"compiler\": " << std::quoted(__VERSION__) << ",\n"
        << "  \"real_digits\": " << std::numeric_limits<Real>::digits << ",\n"
        << "  \"accuracy\": \"LL_plus_single_Oalpha_boundary_FF_only_matching\",\n"
        << "  \"NLL_kernel_included\": false,\n"
        << "  \"hard_Oalpha_included\": false,\n"
        << "  \"running_alpha_included\": false,\n"
        << "  \"mixing_included\": false,\n"
        << "  \"boundary_insertions_max_across_both_legs\": 1,\n"
        << "  \"scheme\": \"MSbar\",\n"
        << "  \"units\": \"GeV_and_dimensionless_acceptances\",\n"
        << "  \"parent_mass\": " << input.mass << ",\n"
        << "  \"muon_mass\": " << input.qed.muon_mass << ",\n"
        << "  \"mu_0\": " << input.qed.mu_0 << ",\n"
        << "  \"mu_frag\": " << input.qed.mu_frag << ",\n"
        << "  \"alpha\": " << input.qed.alpha << ",\n"
        << "  \"time\": " << baseline.evolution().time() << ",\n"
        << "  \"pt_min\": " << input.cuts.pt_min << ",\n"
        << "  \"eta_max\": " << input.cuts.eta_max << ",\n"
        << "  \"mass_min\": " << input.cuts.mass_min << ",\n"
        << "  \"mass_max\": " << input.cuts.mass_max << ",\n"
        << "  \"absolute_tolerance\": " << input.integration.absolute_tolerance << ",\n"
        << "  \"relative_tolerance\": " << input.integration.relative_tolerance << ",\n"
        << "  \"max_intervals\": " << input.integration.max_intervals << ",\n"
        << "  \"born_acceptance\": " << born << ",\n"
        << "  \"fixed_order_FF_acceptance\": " << born+input.qed.expansion_parameter()*fixed << ",\n"
        << "  \"solutions\": {\n    \"baseline\": ";
    write_solution(baseline,input);
    std::cout << ",\n    \"refined_grid\": ";
    write_solution(refined_grid,input);
    std::cout << ",\n    \"refined_cutoff\": ";
    write_solution(refined_cut,input);
    std::cout << "\n  }\n}\n";
  } catch (const std::exception& error) {
    std::cerr << "leptonff_match: " << error.what() << '\n';
    return 1;
  }
}
