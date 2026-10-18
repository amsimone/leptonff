#include "leptonff/decay.hpp"
#include "reference_config.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace leptonff;

int main(int argc, char** argv) {
  try {
    if (argc == 2 && std::string(argv[1]) == "--help") {
      std::cout << "Usage: leptonff_decay [reference.conf]\n"
                   "Outputs a dimensionless isotropic-decay FF reference as JSON.\n";
      return 0;
    }
    if (argc > 2) throw std::invalid_argument("Usage: leptonff_decay [reference.conf]");
    ReferenceInput input;
    if (argc == 2) read_reference_config(argv[1], input);
    const Real mass = input.mass;
    const auto& qed = input.qed;
    const auto& cuts = input.cuts;
    const auto& options = input.integration;
    qed.validate();
    cuts.validate();
    const Real born = isotropic_acceptance(mass, 1, 1, cuts);
    const Real coefficient = decay_ff_coefficient_order1(mass, cuts, qed, options);
    auto tighter = options;
    tighter.absolute_tolerance /= 10;
    tighter.relative_tolerance /= 10;
    const Real refined = decay_ff_coefficient_order1(mass, cuts, qed, tighter);
    const Real a = qed.expansion_parameter();
    const Real rate = born + a*coefficient;
    Real lepton_momentum = 0, photon_momentum = 0;
    for (const auto observed : one_family_species) {
      lepton_momentum += fragmentation_coefficient_order1(
          Species::muon_minus, observed, qed).moment(2, options);
      photon_momentum += fragmentation_coefficient_order1(
          Species::photon, observed, qed).moment(2, options);
    }
    std::cout << std::setprecision(18)
        << "{\n"
        << "  \"schema_version\": 1,\n"
        << "  \"library_version\": \"" << LEPTONFF_VERSION << "\",\n"
        << "  \"build_type\": \"" << LEPTONFF_BUILD_TYPE << "\",\n"
        << "  \"model\": \"isotropic_massless_decay_at_rest\",\n"
        << "  \"accuracy\": \"Born_plus_Oalpha_FF_only\",\n"
        << "  \"hard_Oalpha_included\": false,\n"
        << "  \"resummed\": false,\n"
        << "  \"scheme\": \"MSbar\",\n"
        << "  \"coupling\": \"fixed_alpha\",\n"
        << "  \"families\": 1,\n"
        << "  \"units\": \"GeV_and_dimensionless_acceptances\",\n"
        << "  \"compiler\": " << std::quoted(__VERSION__) << ",\n"
        << "  \"real_digits\": " << std::numeric_limits<Real>::digits << ",\n"
        << "  \"parent_mass\": " << mass << ",\n"
        << "  \"muon_mass\": " << qed.muon_mass << ",\n"
        << "  \"mu_0\": " << qed.mu_0 << ",\n"
        << "  \"mu_frag\": " << qed.mu_frag << ",\n"
        << "  \"alpha\": " << qed.alpha << ",\n"
        << "  \"pt_min\": " << cuts.pt_min << ",\n"
        << "  \"eta_max\": " << cuts.eta_max << ",\n"
        << "  \"mass_min\": " << cuts.mass_min << ",\n"
        << "  \"mass_max\": " << cuts.mass_max << ",\n"
        << "  \"absolute_tolerance\": " << options.absolute_tolerance << ",\n"
        << "  \"relative_tolerance\": " << options.relative_tolerance << ",\n"
        << "  \"max_intervals\": " << options.max_intervals << ",\n"
        << "  \"born_acceptance\": " << born << ",\n"
        << "  \"ff_coefficient_of_alpha_over_2pi\": " << coefficient << ",\n"
        << "  \"born_plus_ff_order1\": " << rate << ",\n"
        << "  \"rate_refinement_shift\": " << a*(refined-coefficient) << ",\n"
        << "  \"lepton_momentum_coefficient_residual\": " << lepton_momentum << ",\n"
        << "  \"photon_momentum_coefficient_residual\": " << photon_momentum << "\n"
        << "}\n";
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "leptonff_decay: " << error.what() << '\n';
    return 1;
  }
}
