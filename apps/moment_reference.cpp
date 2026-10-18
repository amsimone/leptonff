#include "leptonff/moment_evolution.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void print_matrix(const leptonff::MomentMatrix& matrix) {
  std::cout << '[';
  for (unsigned i = 0; i < 3; ++i) {
    if (i) std::cout << ',';
    std::cout << '[';
    for (unsigned j = 0; j < 3; ++j) {
      if (j) std::cout << ',';
      std::cout << matrix[i][j];
    }
    std::cout << ']';
  }
  std::cout << ']';
}
}
int main(int argc, char**) {
  try {
    if (argc != 1) throw std::invalid_argument("Usage: leptonff_moments");
    using namespace leptonff;
    MomentEvolutionParameters p;
    MomentEvolutionOptions controls;
    controls.kernel_integration = {1e-11L,1e-11L,8192};
    std::cout << std::setprecision(18)
      << "{\n\"version\":\"" << LEPTONFF_VERSION << "\",\n"
      << "\"compiler\":\"" << __VERSION__ << "\",\n"
      << "\"real_mantissa_bits\":" << std::numeric_limits<Real>::digits << ",\n"
      << "\"approximation\":\"one-family MSbar real moments; NLO time-like matrix and O(alpha) boundary; inclusive pairs; no x-space inversion\",\n"
      << "\"alpha_msbar_at_mu0\":" << p.alpha_msbar_at_mu0 << ",\n"
      << "\"mu_0\":" << p.mu_0 << ",\n\"mu_frag\":" << p.mu_frag
      << ",\n\"muon_mass\":" << p.muon_mass << ",\n"
      << "\"species\":[\"mu-\",\"mu+\",\"gamma\"],\n"
      << "\"beta0_in_a\":" << 2.0L/3 << ",\n\"beta1_in_a\":1,\n"
      << "\"step_tolerance\":" << controls.step_tolerance << ",\n"
      << "\"kernel_tolerance\":1e-11,\n\"moments\":[\n";
    bool first = true;
    for (auto mode : {QEDCouplingMode::fixed,QEDCouplingMode::one_loop,QEDCouplingMode::two_loop})
      for (Real N : {2.L,3.L,5.L}) {
        p.coupling = mode;
        const auto r = evolve_qed_moment(N,p,controls);
        if (!first) std::cout << ",\n";
        first = false;
        std::cout << "{\"N\":" << N << ",\"coupling\":\"" << coupling_mode_name(mode)
          << "\",\"a_final\":" << r.a_final << ",\"integral_a\":" << r.integral_a
          << ",\"integral_a2\":" << r.integral_a2 << ",\"steps\":" << r.steps
          << ",\"step_refinement_shift\":" << r.step_refinement_shift;
        for (const auto& entry : {std::make_pair("kernel0",r.kernel0),
                                 std::make_pair("kernel1",r.kernel1),
                                 std::make_pair("boundary1",r.boundary1),
                                 std::make_pair("evolution",r.evolution),
                                 std::make_pair("fragmentation",r.fragmentation),
                                 std::make_pair("coefficient1",r.coefficient1),
                                 std::make_pair("coefficient2",r.coefficient2)}) {
          std::cout << ",\"" << entry.first << "\":";
          print_matrix(entry.second);
        }
        std::cout << '}';
      }
    std::cout << "\n]}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
