#include "leptonff/nlo_qed.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

int main(int argc, char**) {
  try {
    if (argc != 1) throw std::invalid_argument("Usage: leptonff_kernels");
    using namespace leptonff;
    using C = NLOLeptonComponent;
    const IntegrationOptions refined{1e-11L,1e-11L,8192};
    Real max_refinement = 0;
    bool first = true;
    std::cout << std::setprecision(18)
      << "{\n  \"version\": \"" << LEPTONFF_VERSION << "\",\n"
      << "  \"compiler\": \"" << __VERSION__ << "\",\n"
      << "  \"real_mantissa_bits\": " << std::numeric_limits<Real>::digits << ",\n"
      << "  \"approximation\": \"one-family MSbar time-like NLO matrix; inclusive pairs; no evolution\",\n"
      << "  \"normalization\": \"coefficient of (alpha/(2*pi))^2\",\n"
      << "  \"photon_channels_available\": true,\n  \"kernels\": [\n";
    auto row = [&](const std::string& name, const Distribution& kernel, bool finite_n1) {
      if (!first) std::cout << ",\n";
      first = false;
      std::cout << "    {\"name\": \"" << name << "\", \"moments\": [";
      for (unsigned N : {1u,2u,3u,5u}) {
        if (N != 1) std::cout << ", ";
        if (N == 1 && !finite_n1) { std::cout << "null"; continue; }
        const Real value = kernel.moment(N,refined);
        max_refinement = std::max(max_refinement,std::abs(value-kernel.moment(N)));
        std::cout << value;
      }
      std::cout << "], \"z2_above_0p37\": "
        << kernel.apply([](Real z) { return z*z; },.37L,refined) << "}";
    };
    for (C component : {C::photonic,C::nonsinglet_pair,C::singlet_pair,
                        C::same_charge_interference,C::opposite_charge_interference})
      row(nlo_lepton_component_name(component),nlo_lepton_component(component),
          component != C::singlet_pair);
    row("same_charge",nlo_lepton_kernel(Species::muon_minus,Species::muon_minus),false);
    row("opposite_charge",nlo_lepton_kernel(Species::muon_minus,Species::muon_plus),false);
    row("valence",nlo_valence_kernel(),true);
    row("lepton_to_photon",next_to_leading_kernel(Species::muon_minus,Species::photon),true);
    row("photon_to_lepton",next_to_leading_kernel(Species::photon,Species::muon_minus),true);
    row("photon_to_photon",next_to_leading_kernel(Species::photon,Species::photon),false);
    Real momentum_residual = 0;
    for (Species parent : one_family_species) {
      Real sum = 0;
      for (Species observed : one_family_species)
        sum += next_to_leading_kernel(parent,observed).moment(2,refined);
      momentum_residual = std::max(momentum_residual,std::abs(sum));
    }
    std::cout << "\n  ],\n  \"moment_indices\": [1,2,3,5],\n"
      << "  \"maximum_moment_refinement_shift\": " << max_refinement << ",\n"
      << "  \"maximum_parent_momentum_residual\": " << momentum_residual << "\n}\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
