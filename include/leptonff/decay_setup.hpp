#pragma once
#include "leptonff/production_decay.hpp"
#include "leptonff/theory.hpp"
#include <string>

namespace leptonff {
enum class DecayPreset {born,photonic_ll,photonic_nll,leptonic_nll,qed_hvp};
DecayPreset decay_preset(const std::string&);
// Presets resolve physics once. Numerical and sampling controls may then be
// set on the returned configuration. Photonic alpha is explicitly fixed;
// the HVP preset requires its on-shell anchor and keeps QCD switched off.
DecayConfiguration make_decay_configuration(SpinDecay,DecayPreset,const AlphaInput&,
    Real muon_mass=.1056583755L,Real final_scale=91.1876L,Real hadronic_scale=2);
// Complete, versioned identity for prepared physics / response persistence.
// Includes numerical controls, thresholds and retained perturbative convention.
std::string decay_configuration_identity(const DecayConfiguration&);
std::string decay_configuration_json(const DecayConfiguration&);
} // namespace leptonff
