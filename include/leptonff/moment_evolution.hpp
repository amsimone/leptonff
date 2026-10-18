#pragma once

#include "leptonff/nlo_qed.hpp"
#include "leptonff/flavours.hpp"
#include <optional>

namespace leptonff {
using MomentMatrix = std::array<std::array<Real,3>,3>;
MomentMatrix identity_moment_matrix();
MomentMatrix multiply_moment_matrices(const MomentMatrix&, const MomentMatrix&);

enum class QEDCouplingMode { fixed, one_loop, two_loop };
const char* coupling_mode_name(QEDCouplingMode mode);

struct QCDEvolutionParameters {
  Real alpha_s_reference=.1193L,reference_scale=91.1876L,start_scale=2;
};
struct HadronicMatching {
  Real scale=2,coefficient=0; // a_above=a_below/(1-a_below*coefficient)
};

struct MomentEvolutionParameters {
  // Explicit MSbar input at mu_0. No automatic on-shell conversion, flavours
  // or thresholds; this value must not silently replace a decay's fixed alpha.
  Real alpha_msbar_at_mu0 = 1/137.035999084L;
  Real muon_mass = .1056583755L;
  Real mu_0 = .1056583755L;
  Real mu_frag = 91.1876L;
  QEDCouplingMode coupling = QEDCouplingMode::two_loop;
  // Opt-in variable-flavour time-like evolution; empty preserves the legacy
  // one-family backend. Output matrices are its mu-/mu+/photon SUBBLOCK,
  // not a closed momentum-conserving three-species system in this mode.
  std::optional<QEDFlavourScheme> flavour_scheme;
  std::optional<QCDEvolutionParameters> qcd;
  std::optional<HadronicMatching> hadronic_matching;
  void validate() const;
};

struct MomentEvolutionOptions {
  Real step_tolerance = 1e-11L; // absolute max-entry shift on global doubling
  std::size_t initial_steps = 16;
  std::size_t max_steps = 65536;
  IntegrationOptions kernel_integration{};
  void validate() const;
};

struct QEDMomentEvolution {
  Real N = 0, log_scale_ratio = 0;
  Real a_initial = 0, a_final = 0; // a=alpha_MSbar/(2*pi)
  Real integral_a = 0, integral_a2 = 0;
  MomentMatrix kernel0{}, kernel1{}, boundary1{};
  MomentMatrix evolution{}, fragmentation{};
  // D = I + a_initial*coefficient1 + a_initial^2*coefficient2 + ... .
  // Unknown constant two-loop boundary omitted; retained logarithms explicit.
  MomentMatrix coefficient1{}, coefficient2{};
  // All parent rows, selected observed columns {mu-,mu+,gamma}; opt-in only.
  std::vector<std::array<Real,3>> flavour_fragmentation,flavour_evolution,
      flavour_coefficient1,flavour_coefficient2;
  std::size_t steps = 0;
  Real step_refinement_shift = 0;
};

// Real N>1 coupled moment reference. Solves dE/dL=(a P0+a^2 P1) E,
// da/dL=2a^2/3+a^3 in two-loop mode; D=E (I+a_initial*d1).
// Step-doubled RK4 keeps matrix ordering. This is not an x-space inversion,
// a fiducial prediction or an exclusive extra-pair prescription.
QEDMomentEvolution evolve_qed_moment(Real N,
    const MomentEvolutionParameters& parameters = {},
    const MomentEvolutionOptions& options = {});
struct FlavourRunningResult {
  Real a_final=0,integral_a=0,integral_a2=0,beta0_integral=0;
};
// Piecewise two-loop MSbar flow and one-loop decoupling, forward or backward.
// beta0_integral includes threshold logs: a(Q)=a0+a0^2*B(Q)+... .
FlavourRunningResult run_flavour_coupling(const MomentEvolutionParameters&,Real scale);
// LO QCD running, continuous across heavy thresholds at this accuracy.
Real run_qcd_coupling(const MomentEvolutionParameters&,Real scale);
Real integral_qcd_coupling(const MomentEvolutionParameters&,Real lower,Real upper);
} // namespace leptonff
