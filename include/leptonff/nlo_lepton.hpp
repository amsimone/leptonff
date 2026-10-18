#pragma once

#include "leptonff/fragmentation.hpp"

namespace leptonff {
// One-family, unit-charge, MSbar time-like coefficients of (alpha/(2*pi))^2.
// These are inclusive identified-lepton densities, including pair production.
// Source-separated terms are not independently number conserving.
enum class NLOLeptonComponent {
  photonic,
  nonsinglet_pair,
  singlet_pair,
  same_charge_interference,
  opposite_charge_interference
};
const char* nlo_lepton_component_name(NLOLeptonComponent component);
Distribution nlo_lepton_component(NLOLeptonComponent component);

// Only the charged-lepton block is available. Photon parents/daughters throw;
// a missing entry must never masquerade as a zero kernel. Not an NLL solver.
Distribution nlo_lepton_kernel(Species parent, Species observed);

// Same minus opposite charge, with the singular singlet cancelled analytically
// BEFORE constructing a Distribution. Its N=1 charge moment is finite and zero.
Distribution nlo_valence_kernel();
} // namespace leptonff
