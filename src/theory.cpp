#include "leptonff/theory.hpp"
#include <cmath>
#include <stdexcept>

namespace leptonff {
AlphaInput AlphaInput::on_shell(Real value) { return {value, AlphaScheme::on_shell, 0}; }
AlphaInput AlphaInput::msbar(Real value, Real scale) { return {value, AlphaScheme::msbar, scale}; }

void AlphaInput::validate() const {
  if (!std::isfinite(value) || value < 0 || value >= 1 ||
      !std::isfinite(reference_scale))
    throw std::invalid_argument("Invalid alpha input");
  if (scheme == AlphaScheme::on_shell) {
    if (reference_scale != 0)
      throw std::invalid_argument("On-shell alpha has no reference scale; use AlphaInput::on_shell");
  } else if (scheme == AlphaScheme::msbar) {
    if (reference_scale <= 0)
      throw std::invalid_argument("MSbar alpha needs a positive reference scale");
  } else throw std::invalid_argument("Unknown alpha scheme");
}

MomentEvolutionParameters EvolutionTheory::parameters() const {
  alpha.validate();
  MomentEvolutionParameters p;
  p.muon_mass=muon_mass; p.mu_0=initial_scale; p.mu_frag=final_scale;
  p.coupling=running; p.flavour_scheme=flavours;
  if (alpha.scheme == AlphaScheme::msbar) {
    if (alpha.reference_scale != initial_scale)
      throw std::invalid_argument("MSbar alpha must be provided at the initial scale; run it explicitly first");
    p.alpha_msbar_at_mu0=alpha.value;
  } else {
    const auto scheme=flavours.value_or(QEDFlavourScheme{{{"muon",-1,muon_mass,1,0,false}},0,false});
    scheme.validate();
    p.alpha_msbar_at_mu0=leptonic_alpha_msbar(alpha.value,initial_scale,scheme);
  }
  p.validate();
  return p;
}
} // namespace leptonff
