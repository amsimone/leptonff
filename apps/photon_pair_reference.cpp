#include "leptonff/photon_pair.hpp"
#include "reference_config.hpp"
#include <iomanip>
#include <iostream>
#include <limits>

using namespace leptonff;
namespace {
void coefficient(const PhotonPairHardCoefficient& r) {
  std::cout << "{\"projected\":" << r.projected << ",\"geometry\":" << r.geometry
      << ",\"total\":" << r.total << '}';
}
}
int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("Usage: leptonff_photon_pair config");
    ReferenceInput input;
    read_reference_config(argv[1],input);
    const Real M = input.mass;
    HardDecayOptions baseline;
    baseline.integration = {1e-7L,1e-7L,8192};
    baseline.collinear_cutoff = 1e-6L;
    auto refined = baseline;
    refined.collinear_cutoff = 1e-8L;
    refined.integration = {2.5e-8L,2.5e-8L,8192};
    auto final = refined;
    final.collinear_cutoff = 2.5e-9L;
    final.integration = {6.25e-9L,6.25e-9L,16384};
    std::cout << std::setprecision(18)
        << "{\"schema_version\":1,\"library_version\":\"" << LEPTONFF_VERSION
        << "\",\"compiler\":" << std::quoted(__VERSION__)
        << ",\"real_digits\":" << std::numeric_limits<Real>::digits
        << ",\"model\":\"unpolarised_vector_at_rest_one_specified_primary_charge\""
        << ",\"normalisation\":\"coefficient_of_alpha_over_2pi_relative_to_massless_vector_Born\""
        << ",\"scheme\":\"MSbar\",\"evolved_FF_convolution_included\":false"
        << ",\"parent_mass\":" << M << ",\"muon_mass\":" << input.qed.muon_mass
        << ",\"mu_frag\":" << input.qed.mu_frag
        << ",\"pt_min\":" << input.cuts.pt_min << ",\"eta_max\":" << input.cuts.eta_max
        << ",\"mass_min\":" << input.cuts.mass_min << ",\"mass_max\":" << input.cuts.mass_max
        << ",\"collinear_cutoff_final\":" << final.collinear_cutoff
        << ",\"density_scan\":[";
    bool first = true;
    for (Real y : {.2L,.6L,.85L}) for (Real factor : {1.0L,.5L,.25L}) {
      auto q = input.qed;
      q.muon_mass *= factor; q.mu_0 = q.muon_mass;
      const Real H = vector_photon_pair_mass_hard_coefficient(M,q.mu_frag).resolved_density(y);
      const Real ff = fragmentation_coefficient_order1(Species::muon_minus,
          Species::photon,q).resolved_density(y);
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"y\":" << y << ",\"muon_mass\":" << q.muon_mass
          << ",\"hard\":" << H << ",\"ff\":" << ff << ",\"leading_power\":" << H+ff
          << ",\"massive\":" << massive_vector_photon_pair_spectrum(M,q.muon_mass,y) << '}';
    }
    std::cout << "],\"moments\":[";
    first = true;
    for (Real N : {2.0L,3.0L,5.0L}) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"N\":" << N << ",\"hard\":"
          << vector_photon_pair_mass_hard_coefficient(M,input.qed.mu_frag).moment(N)
          << '}';
    }
    std::cout << "],\"fiducial\":[";
    first = true;
    for (auto fractions : {std::pair<Real,Real>{1,1},{.9L,.95L},{.75L,.8L}}) {
      if (!first) std::cout << ',';
      first = false;
      const auto [fg,fm] = fractions;
      std::cout << "{\"photon_fraction\":" << fg << ",\"muon_fraction\":" << fm
          << ",\"baseline\":";
      coefficient(fiducial_vector_photon_pair_hard(M,input.qed.mu_frag,input.cuts,fg,fm,baseline));
      std::cout << ",\"refined\":";
      coefficient(fiducial_vector_photon_pair_hard(M,input.qed.mu_frag,input.cuts,fg,fm,refined));
      std::cout << ",\"final\":";
      coefficient(fiducial_vector_photon_pair_hard(M,input.qed.mu_frag,input.cuts,fg,fm,final));
      std::cout << '}';
    }
    std::cout << "]}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n'; return 1;
  }
}
