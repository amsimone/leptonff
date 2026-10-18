#include "leptonff/fiducial_vector.hpp"
#include "leptonff/evolution.hpp"
#include "reference_config.hpp"
#include "ll_grid_options.hpp"

#include <iomanip>
#include <iostream>
#include <limits>

using namespace leptonff;
namespace {
void write_fixed(const FiducialVectorOrder1& r) {
  std::cout << "{\"born\":" << r.born
      << ",\"projected_hard_coefficient\":" << r.projected_hard_coefficient
      << ",\"geometry_coefficient\":" << r.geometry_coefficient
      << ",\"hard_coefficient\":" << r.hard_coefficient
      << ",\"ff_coefficient\":" << r.ff_coefficient
      << ",\"full_coefficient\":" << r.full_coefficient
      << ",\"rate\":" << r.rate << '}';
}
void write_massive(const MassiveFiducialVector& r) {
  std::cout << "{\"born\":" << r.born << ",\"coefficient\":" << r.coefficient
      << ",\"rate\":" << r.rate << '}';
}
void write_controls(const HardDecayOptions& o) {
  std::cout << "{\"soft_cutoff\":" << o.soft_cutoff
      << ",\"collinear_cutoff\":" << o.collinear_cutoff
      << ",\"absolute_tolerance\":" << o.integration.absolute_tolerance
      << ",\"relative_tolerance\":" << o.integration.relative_tolerance
      << ",\"max_intervals\":" << o.integration.max_intervals << '}';
}
void write_solution(const ReferenceInput& input, const LLGrid& grid,
                    const FiducialVectorOrder1& fixed) {
  const LLBoundaryNonSinglet evolution(input.qed,grid);
  const auto result = fiducial_vector_ll_match(input.mass,input.cuts,evolution,fixed,input.integration);
  const auto& d = evolution.evolution().diagnostics();
  std::cout << "{\"bins\":" << grid.bins
      << ",\"log_span\":" << grid.log_span
      << ",\"spacing\":" << grid.spacing()
      << ",\"soft_log_cutoff\":" << grid.soft_log_cutoff
      << ",\"wrap_probability_bound\":" << d.wrap_probability_bound
      << ",\"ff_matched_rate\":" << result.ff_matched
      << ",\"hard_addition\":" << result.hard_addition
      << ",\"matched_rate\":" << result.matched
      << ",\"fixed_order_rate\":" << result.fixed_order << '}';
}
} // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << "Usage: leptonff_fiducial config [bins [soft_log_cutoff [log_span]]]\n"
                   "Unpolarised vector at rest: full Oalpha plus additive LL/boundary matching.\n";
      return 0;
    }
    if (argc < 2 || argc > 5)
      throw std::invalid_argument("Usage: leptonff_fiducial config [bins [soft_log_cutoff [log_span]]]");
    ReferenceInput input;
    read_reference_config(argv[1],input);
    const auto grid = read_ll_grid(argc,argv);
    HardDecayOptions baseline_options;
    const auto baseline = fiducial_vector_order1(input.mass,input.cuts,input.qed,
                                                baseline_options,input.integration);
    auto soft_options = baseline_options;
    soft_options.soft_cutoff /= 2;
    const auto soft = fiducial_vector_order1(input.mass,input.cuts,input.qed,soft_options,input.integration);
    auto collinear_options = soft_options;
    collinear_options.collinear_cutoff /= 2;
    const auto collinear = fiducial_vector_order1(input.mass,input.cuts,input.qed,
                                                 collinear_options,input.integration);
    auto final_options = collinear_options;
    final_options.integration.absolute_tolerance /= 4;
    final_options.integration.relative_tolerance /= 4;
    const auto final = fiducial_vector_order1(input.mass,input.cuts,input.qed,
                                             final_options,input.integration);
    const auto massive_baseline = massive_fiducial_vector_order1(input.mass,input.cuts,
                                                                 input.qed,baseline_options);
    const auto massive = massive_fiducial_vector_order1(input.mass,input.cuts,input.qed,final_options);
    std::cout << std::setprecision(18)
        << "{\n  \"schema_version\":1,\n  \"library_version\":\"" << LEPTONFF_VERSION
        << "\",\n  \"build_type\":\"" << LEPTONFF_BUILD_TYPE
        << "\",\n  \"compiler\":" << std::quoted(__VERSION__)
        << ",\n  \"real_digits\":" << std::numeric_limits<Real>::digits
        << ",\n  \"model\":\"unpolarised_vector_at_rest_with_symmetric_fiducial_cuts\",\n"
        << "  \"accuracy\":\"Oalpha_plus_LL_and_single_evolved_boundary_additive_matching\",\n"
        << "  \"normalisation\":\"massless_vector_Born\",\n"
        << "  \"scheme\":\"MSbar\",\n  \"coupling\":\"fixed_alpha\",\n"
        << "  \"hard_Oalpha_included\":true,\n  \"hard_LL_dressing_included\":false,\n"
        << "  \"NLL_kernel_included\":false,\n  \"running_alpha_included\":false,\n"
        << "  \"mixing_included\":false,\n  \"axial_current\":false,\n"
        << "  \"boundary_insertions_max_across_both_legs\":1,\n"
        << "  \"parent_mass\":" << input.mass
        << ",\n  \"muon_mass\":" << input.qed.muon_mass
        << ",\n  \"mu_0\":" << input.qed.mu_0
        << ",\n  \"mu_frag\":" << input.qed.mu_frag
        << ",\n  \"alpha\":" << input.qed.alpha
        << ",\n  \"pt_min\":" << input.cuts.pt_min
        << ",\n  \"eta_max\":" << input.cuts.eta_max
        << ",\n  \"mass_min\":" << input.cuts.mass_min
        << ",\n  \"mass_max\":" << input.cuts.mass_max
        << ",\n  \"absolute_tolerance\":" << input.integration.absolute_tolerance
        << ",\n  \"relative_tolerance\":" << input.integration.relative_tolerance
        << ",\n  \"hard_controls_baseline\":";
    write_controls(baseline_options);
    std::cout << ",\n  \"hard_controls_final\":";
    write_controls(final_options);
    std::cout << ",\n  \"fixed_baseline\":";
    write_fixed(baseline);
    std::cout << ",\n  \"fixed\":";
    write_fixed(final);
    std::cout << ",\n  \"hard_soft_refinement_rate_shift\":" << soft.rate-baseline.rate
        << ",\n  \"hard_collinear_refinement_rate_shift\":" << collinear.rate-soft.rate
        << ",\n  \"hard_tolerance_refinement_rate_shift\":" << final.rate-collinear.rate
        << ",\n  \"massive\":";
    write_massive(massive);
    std::cout << ",\n  \"massive_refinement_rate_shift\":" << massive.rate-massive_baseline.rate
        << ",\n  \"scale_scan\":[";
    bool first = true;
    for (Real factor : {.5L,1.0L,2.0L}) {
      auto q = input.qed;
      q.mu_frag = input.mass*factor;
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"mu_frag\":" << q.mu_frag << ",\"fixed\":";
      write_fixed(fiducial_vector_order1(input.mass,input.cuts,q,final_options,input.integration));
      std::cout << '}';
    }
    std::cout << "],\n  \"mass_scan\":[";
    first = true;
    for (Real factor : {1.0L,.5L,.25L}) {
      auto q = input.qed;
      q.muon_mass *= factor;
      q.mu_0 *= factor;
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"muon_mass\":" << q.muon_mass << ",\"fixed\":";
      write_fixed(fiducial_vector_order1(input.mass,input.cuts,q,final_options,input.integration));
      std::cout << ",\"massive\":";
      write_massive(massive_fiducial_vector_order1(input.mass,input.cuts,q,final_options));
      std::cout << '}';
    }
    std::cout << "],\n  \"solutions\":{\"baseline\":";
    write_solution(input,grid,final);
    auto fine_grid = grid;
    fine_grid.bins *= 2;
    std::cout << ",\"refined_grid\":";
    write_solution(input,fine_grid,final);
    fine_grid.soft_log_cutoff /= 2;
    std::cout << ",\"refined_cutoff\":";
    write_solution(input,fine_grid,final);
    std::cout << "}\n}\n";
  } catch (const std::exception& e) {
    std::cerr << "leptonff_fiducial: " << e.what() << '\n';
    return 1;
  }
}
