#pragma once
#include "leptonff/moment_evolution.hpp"

namespace leptonff {
enum class AlphaScheme { on_shell, msbar };

// Running order is separate: QEDCouplingMode belongs to EvolutionTheory.
// An on-shell input has reference_scale=0; an MSbar input names its scale.
struct AlphaInput {
  Real value = 1/137.035999084L;
  AlphaScheme scheme = AlphaScheme::msbar;
  Real reference_scale = .1056583755L;
  static AlphaInput on_shell(Real value);
  static AlphaInput msbar(Real value, Real scale);
  void validate() const;
};

// Physical settings only. Grids and solver tolerances remain in LLGrid,
// QEDXSpaceGrid and MomentEvolutionOptions. The legacy explicit backend
// parameters remain available for advanced QCD/HVP studies.
struct EvolutionTheory {
  AlphaInput alpha;
  Real muon_mass = .1056583755L;
  Real initial_scale = .1056583755L;
  Real final_scale = 91.1876L;
  QEDCouplingMode running = QEDCouplingMode::two_loop;
  std::optional<QEDFlavourScheme> flavours;
  MomentEvolutionParameters parameters() const;
};
} // namespace leptonff
