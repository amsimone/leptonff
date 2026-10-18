#pragma once

#include "leptonff/distribution.hpp"

#include <array>

namespace leptonff {
enum class Species { muon_minus, muon_plus, photon };
constexpr std::array<Species, 3> one_family_species{
    Species::muon_minus, Species::muon_plus, Species::photon};
const char* species_name(Species species);

struct QEDParameters {
  Real alpha = 1 / 137.035999084L;
  Real muon_mass = .1056583755L; // GeV
  Real mu_0 = .1056583755L;     // GeV
  Real mu_frag = 91.1876L;     // GeV
  void validate() const;
  Real expansion_parameter() const;
};

// Coefficient of alpha/(2*pi), NOT the full FF. MSbar, one family, unit charge.
// Frixione 1909.03886v2: eqs. 289,293,300,302 -> 142,210,255,260.
Distribution initial_condition_coefficient(Species parent, Species observed,
                                            const QEDParameters& parameters);

// Time-like LO kernels in the explicit parent -> observed convention.
Distribution leading_kernel(Species parent, Species observed);

// D1(mu_frag) = d1(mu_0) + log(mu_frag^2/mu_0^2) P0.
// Fixed-order scale transport only; no all-order evolution or NLL claim.
Distribution fragmentation_coefficient_order1(Species parent, Species observed,
                                              const QEDParameters& parameters);
Distribution fragmentation_order1(Species parent, Species observed,
                                   const QEDParameters& parameters);
} // namespace leptonff
