#include "leptonff/hard_dressing.hpp"
#include "reference_config.hpp"
#include "ll_grid_options.hpp"

#include <iomanip>
#include <iostream>
#include <limits>

using namespace leptonff;
namespace {
void head(const std::vector<Real>& v) {
  std::cout << '[';
  for (std::size_t j=0;j<std::min(std::size_t{16},v.size());++j) {
    if (j) std::cout << ',';
    std::cout << v[j];
  }
  std::cout << ']';
}
void solution(const ReferenceInput& in, const LLGrid& grid,
              const HardDressingOptions& o, const FiducialVectorOrder1& fixed) {
  const LLBoundaryNonSinglet boundary(in.qed,grid);
  const auto r=fiducial_vector_hard_dressed_match(in.mass,in.cuts,boundary,fixed,o,in.integration);
  const auto& ll=boundary.evolution(); const auto& p=r.projection; const auto& q=r.quadrature;
  std::cout << "{\"bins\":" << grid.bins << ",\"spacing\":" << grid.spacing()
      << ",\"log_span\":" << grid.log_span << ",\"soft_log_cutoff\":" << grid.soft_log_cutoff
      << ",\"ff_panels\":" << o.ff_panels << ",\"hard_soft_cutoff\":" << o.geometry.soft_cutoff
      << ",\"hard_collinear_cutoff\":" << o.geometry.collinear_cutoff
      << ",\"hard_absolute_tolerance\":" << o.geometry.integration.absolute_tolerance
      << ",\"hard_relative_tolerance\":" << o.geometry.integration.relative_tolerance
      << ",\"hard_max_intervals\":" << o.geometry.integration.max_intervals
      << ",\"acceptance\":" << decay_ll_acceptance(in.mass,in.cuts,ll)
      << ",\"resolved_jump_rate\":" << ll.diagnostics().resolved_jump_rate
      << ",\"LL_wrap_probability_bound\":" << ll.diagnostics().wrap_probability_bound
      << ",\"hard_wrap_coefficient_bound\":" << p.wrap_coefficient_bound
      << ",\"hard_omitted_tail_coefficient_bound\":" << p.omitted_tail_coefficient_bound
      << ",\"hard_number_residual\":" << p.number_residual
      << ",\"hard_imaginary_weight_mass\":" << p.imaginary_weight_mass
      << ",\"quadrature_number_residual\":" << q.number_residual
      << ",\"quadrature_first_moment_residual\":" << q.first_moment_residual
      << ",\"quadrature_negative_weight_mass\":" << q.negative_weight_mass
      << ",\"ff_matched_rate\":" << r.ff_matched
      << ",\"undressed_hard_coefficient\":" << r.undressed_hard_coefficient
      << ",\"projected_dressing_shift\":" << r.projected_dressing_shift
      << ",\"geometry_dressing_shift\":" << r.geometry_dressing_shift
      << ",\"dressed_hard_coefficient\":" << r.dressed_hard_coefficient
      << ",\"additive_matched_rate\":" << r.additive_matched
      << ",\"matched_rate\":" << r.matched << ",\"weights_head\":";
  head(ll.weights()); std::cout << ",\"initial_hard_head\":";
  head(p.initial_weights); std::cout << ",\"evolved_hard_head\":";
  head(p.evolved_weights); std::cout << ",\"quadrature_nodes\":[";
  bool first=true;
  for (const auto& node:q.nodes) {
    if (!first) std::cout << ','; first=false;
    std::cout << "{\"fraction\":" << node.fraction << ",\"weight\":" << node.weight << '}';
  }
  std::cout << "]}";
}
} // namespace

int main(int argc,char** argv) {
  try {
    if (argc==2 && std::string(argv[1])=="--help") {
      std::cout << "Usage: leptonff_dressed config [bins [LL_cutoff [log_span [FF_panels]]]]\n"
                   "Dresses H1 with pure LL on both muons; one finite insertion globally.\n";
      return 0;
    }
    if (argc<2 || argc>6) throw std::invalid_argument("Usage: leptonff_dressed config [bins [LL_cutoff [log_span [FF_panels]]]]");
    ReferenceInput in; read_reference_config(argv[1],in);
    const LLGrid grid=read_ll_grid(argc,argv);
    HardDressingOptions controls; controls.ff_panels=16;
    if (argc==6) {
      const Real panels=parse_ll_real(argv[5]);
      if (panels<1 || panels>128 || std::floor(panels)!=panels)
        throw std::invalid_argument("Driver FF panels must be integers in [1,128] for refinement");
      controls.ff_panels=static_cast<std::size_t>(panels);
    }
    controls.validate();
    auto fc=controls.geometry; fc.soft_cutoff/=2; fc.collinear_cutoff/=2;
    fc.integration.absolute_tolerance/=4; fc.integration.relative_tolerance/=4;
    const auto fixed=fiducial_vector_order1(in.mass,in.cuts,in.qed,fc,in.integration);
    std::cout << std::setprecision(18)
        << "{\n\"schema_version\":1,\n\"library_version\":\"" << LEPTONFF_VERSION
        << "\",\n\"build_type\":\"" << LEPTONFF_BUILD_TYPE << "\",\n\"compiler\":" << std::quoted(__VERSION__)
        << ",\n\"real_digits\":" << std::numeric_limits<Real>::digits
        << ",\n\"model\":\"unpolarised_vector_at_rest_with_symmetric_fiducial_cuts\",\n"
        << "\"accuracy\":\"LL_plus_single_finite_boundary_or_LL_dressed_H1_with_Oalpha_matching\",\n"
        << "\"normalisation\":\"massless_vector_Born\",\n\"scheme\":\"MSbar\",\n\"coupling\":\"fixed_alpha\",\n"
        << "\"hard_LL_dressing_included\":true,\n\"NLL_kernel_included\":false,\n"
        << "\"hard_boundary_products_included\":false,\n\"running_alpha_included\":false,\n"
        << "\"mixing_included\":false,\n\"axial_current\":false,\n"
        << "\"parent_mass\":" << in.mass << ",\n\"muon_mass\":" << in.qed.muon_mass
        << ",\n\"mu_0\":" << in.qed.mu_0 << ",\n\"mu_frag\":" << in.qed.mu_frag << ",\n\"alpha\":" << in.qed.alpha
        << ",\n\"time\":" << in.qed.expansion_parameter()*2*std::log(in.qed.mu_frag/in.qed.mu_0)
        << ",\n\"pt_min\":" << in.cuts.pt_min << ",\n\"eta_max\":" << in.cuts.eta_max
        << ",\n\"mass_min\":" << in.cuts.mass_min << ",\n\"mass_max\":" << in.cuts.mass_max
        << ",\n\"absolute_tolerance\":" << in.integration.absolute_tolerance
        << ",\n\"relative_tolerance\":" << in.integration.relative_tolerance
        << ",\n\"fixed_order_rate\":" << fixed.rate << ",\n\"fixed_hard_coefficient\":" << fixed.hard_coefficient
        << ",\n\"fixed_hard_soft_cutoff\":" << fc.soft_cutoff
        << ",\n\"fixed_hard_collinear_cutoff\":" << fc.collinear_cutoff
        << ",\n\"fixed_hard_absolute_tolerance\":" << fc.integration.absolute_tolerance
        << ",\n\"fixed_hard_relative_tolerance\":" << fc.integration.relative_tolerance
        << ",\n\"solutions\":{\"baseline\":";
    solution(in,grid,controls,fixed);
    controls.ff_panels*=2; std::cout << ",\n\"refined_ff_quadrature\":";
    solution(in,grid,controls,fixed);
    auto fine=grid; fine.bins*=2; std::cout << ",\n\"refined_grid\":";
    solution(in,fine,controls,fixed);
    fine.soft_log_cutoff/=2; std::cout << ",\n\"refined_cutoff\":";
    solution(in,fine,controls,fixed);
    controls.geometry=fc; std::cout << ",\n\"refined_hard_geometry\":";
    solution(in,fine,controls,fixed);
    std::cout << "}\n}\n";
  } catch (const std::exception& e) {
    std::cerr << "leptonff_dressed: " << e.what() << '\n'; return 1;
  }
}
