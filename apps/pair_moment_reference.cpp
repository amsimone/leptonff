#include "leptonff/pair_moment.hpp"
#include "reference_config.hpp"
#include <iomanip>
#include <iostream>
#include <limits>

int main(int argc, char** argv) {
  try {
    using namespace leptonff;
    if (argc != 2) throw std::invalid_argument("Usage: leptonff_pair_moments config");
    ReferenceInput input;
    read_reference_config(argv[1],input);
    MomentEvolutionParameters p;
    p.alpha_msbar_at_mu0 = input.qed.alpha;
    p.muon_mass = input.qed.muon_mass; p.mu_0 = input.qed.mu_0;
    MomentEvolutionOptions options;
    options.kernel_integration = {1e-10L,1e-10L,8192};
    std::cout << std::setprecision(18)
        << "{\"version\":\"" << LEPTONFF_VERSION << "\",\"compiler\":" << std::quoted(__VERSION__)
        << ",\"real_digits\":" << std::numeric_limits<Real>::digits
        << ",\"model\":\"one_family_unpolarised_vector_inclusive_distinct_parent_pairs\""
        << ",\"accuracy\":\"two_loop_rows_Oalpha_hard_matched_NLL_other_rows_running_diagnostics\""
        << ",\"fiducial_cuts_applied\":false,\"scheme\":\"MSbar\""
        << ",\"parent_mass\":" << input.mass << ",\"mu_0\":" << p.mu_0
        << ",\"muon_mass\":" << p.muon_mass << ",\"alpha_msbar_at_mu0\":" << p.alpha_msbar_at_mu0
        << ",\"hard_coupling_scale\":" << input.mass << ",\"moments\":[";
    bool first = true;
    for (auto mode : {QEDCouplingMode::fixed,QEDCouplingMode::one_loop,QEDCouplingMode::two_loop})
      for (Real N : {2.0L,3.0L,5.0L}) for (Real scale : {.5L,1.0L,2.0L}) {
        p.coupling = mode; p.mu_frag = scale*input.qed.mu_frag;
        const auto r = vector_pair_moment(input.mass,N,p,options);
        if (!first) std::cout << ',';
        first = false;
        std::cout << "{\"N\":" << N << ",\"coupling\":" << std::quoted(coupling_mode_name(mode))
            << ",\"mu_frag\":" << p.mu_frag << ",\"a_hard\":" << r.a_hard
            << ",\"ff\":" << r.ff << ",\"muon_hard\":" << r.muon_hard
            << ",\"photon_hard\":" << r.photon_hard << ",\"matched\":" << r.matched
            << ",\"hard_muon_coefficient\":" << r.hard_muon_coefficient
            << ",\"hard_photon_coefficient\":" << r.hard_photon_coefficient
            << ",\"coefficient1\":" << r.coefficient1 << ",\"coefficient2\":" << r.coefficient2 << '}';
      }
    std::cout << "]}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
