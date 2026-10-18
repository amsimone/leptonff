#include "leptonff/moment_evolution.hpp"
#include "matrix_evolution.hpp"
#include "flavour_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff {
MomentMatrix identity_moment_matrix() {
  MomentMatrix result{};
  for (unsigned i = 0; i < 3; ++i) result[i][i] = 1;
  return result;
}

MomentMatrix multiply_moment_matrices(const MomentMatrix& a, const MomentMatrix& b) {
  MomentMatrix result{};
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      for (unsigned k = 0; k < 3; ++k) result[i][j] += a[i][k]*b[k][j];
  return result;
}

const char* coupling_mode_name(QEDCouplingMode mode) {
  switch (mode) {
  case QEDCouplingMode::fixed: return "fixed";
  case QEDCouplingMode::one_loop: return "one_loop";
  case QEDCouplingMode::two_loop: return "two_loop";
  }
  throw std::invalid_argument("Unknown QED coupling mode");
}

void MomentEvolutionParameters::validate() const {
  QEDParameters{alpha_msbar_at_mu0,muon_mass,mu_0,mu_frag}.validate();
  coupling_mode_name(coupling);
  if(qcd) {
    if(!flavour_scheme || !flavour_scheme->include_gluon || !std::isfinite(qcd->alpha_s_reference) ||
       qcd->alpha_s_reference<0 || qcd->alpha_s_reference>=1 || !std::isfinite(qcd->start_scale) ||
       qcd->start_scale<1 || !std::isfinite(qcd->reference_scale) || qcd->reference_scale<qcd->start_scale)
      throw std::invalid_argument("Invalid LO QCD evolution input");
    for(const auto& f:flavour_scheme->fermions)if(f.colours!=1 && (f.colours!=3 || f.threshold()<qcd->start_scale))
      throw std::invalid_argument("QCD parents require Nc=3 and thresholds at/above the perturbative start");
  }
  if(hadronic_matching && (!flavour_scheme || coupling==QEDCouplingMode::fixed ||
     !std::isfinite(hadronic_matching->scale) || hadronic_matching->scale<1 ||
     !std::isfinite(hadronic_matching->coefficient) || std::abs(hadronic_matching->coefficient)>100))
    throw std::invalid_argument("Invalid hadronic matching input");
  if(flavour_scheme) {
    flavour_scheme->validate();const auto& observed=flavour_scheme->fermions[flavour_scheme->muon_flavour];
    if(std::abs(observed.mass-muon_mass)>1e-12L*muon_mass || observed.threshold()>mu_0)
      throw std::invalid_argument("Observed muon mass/initial active basis is inconsistent with flavour scheme");
  }
}

void MomentEvolutionOptions::validate() const {
  if (!std::isfinite(step_tolerance) || step_tolerance <= 0 ||
      initial_steps == 0 || max_steps > 1048576 || initial_steps > max_steps/2)
    throw std::invalid_argument("Invalid moment evolution refinement controls");
}


QEDMomentEvolution evolve_qed_moment(Real N, const MomentEvolutionParameters& parameters,
                                    const MomentEvolutionOptions& options) {
  parameters.validate();
  options.validate();
  if (!std::isfinite(N) || N <= 1)
    throw std::invalid_argument("Coupled QED moments require finite N>1");
  if(parameters.flavour_scheme)return detail::evolve_flavour_moment(N,parameters,options);
  QEDMomentEvolution result;
  result.N = N;
  result.log_scale_ratio = 2*(std::log(parameters.mu_frag)-std::log(parameters.mu_0));
  result.a_initial = parameters.alpha_msbar_at_mu0/(2*std::acos(-1.0L));
  const QEDParameters boundary_parameters{parameters.alpha_msbar_at_mu0,
      parameters.muon_mass,parameters.mu_0,parameters.mu_frag};
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      const auto parent = one_family_species[i], observed = one_family_species[j];
      result.kernel0[i][j] = leading_kernel(parent,observed).moment(N,options.kernel_integration);
      result.kernel1[i][j] = next_to_leading_kernel(parent,observed).moment(N,options.kernel_integration);
      result.boundary1[i][j] = initial_condition_coefficient(parent,observed,boundary_parameters)
          .moment(N,options.kernel_integration);
    }

  const auto state = detail::evolve_generators(result.kernel0,result.kernel1,
      result.a_initial,result.log_scale_ratio,parameters.coupling,options,
      result.steps,result.step_refinement_shift);
  result.a_final = state.a;
  result.integral_a = state.integral_a;
  result.integral_a2 = state.integral_a2;
  result.evolution = state.e;
  MomentMatrix boundary = identity_moment_matrix();
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      boundary[i][j] += result.a_initial*result.boundary1[i][j];
  result.fragmentation = multiply_moment_matrices(result.evolution,boundary);
  const auto p0_squared = multiply_moment_matrices(result.kernel0,result.kernel0);
  const auto p0_boundary = multiply_moment_matrices(result.kernel0,result.boundary1);
  const Real L = result.log_scale_ratio;
  const Real beta0 = parameters.coupling == QEDCouplingMode::fixed ? 0 : 2.0L/3;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      result.coefficient1[i][j] = result.boundary1[i][j]+L*result.kernel0[i][j];
      result.coefficient2[i][j] = L*result.kernel1[i][j]+L*p0_boundary[i][j]
          +L*L/2*(p0_squared[i][j]+beta0*result.kernel0[i][j]);
    }
  return result;
}
} // namespace leptonff
