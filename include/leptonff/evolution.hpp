#pragma once

#include "leptonff/quadrature.hpp"
#include "leptonff/fragmentation.hpp"

namespace leptonff {
struct LLGrid {
  std::size_t bins = 262144;       // power of two; periodic log-loss lattice
  Real log_span = 32;
  Real soft_log_cutoff = 2e-4L;    // unresolved u=-log(z), not a physics parameter
  Real max_wrap_probability = 1e-10L;
  void validate() const;
  Real spacing() const;
};

struct LLDiagnostics {
  Real resolved_jump_rate = 0;
  Real no_resolved_emission_weight = 1;
  Real wrap_probability_bound = 0;
  Real omitted_large_jump_probability_bound = 0;
  Real mean_log_loss_discretisation_bound = 0;
  Real normalisation_residual = 0;
  Real negative_weight_mass = 0;
  Real imaginary_weight_mass = 0;
};

struct LatticeConvolution {
  std::vector<Real> weights;
  Real imaginary_weight_mass = 0;
  Real number_residual = 0;
};
// Signed periodic convolution; both arrays must have the same power-of-two size.
// The caller is responsible for a physical bound on periodic wrap.
LatticeConvolution circular_lattice_convolution(const std::vector<Real>& first,
                                               const std::vector<Real>& second);

// P0(N) for the whole-plus kernel [(1+z^2)/(1-z)]+; real N>0.
Real non_singlet_kernel_moment(Real N);

// Restricted fixed-alpha photonic/non-singlet LL evolution from delta(1-z).
// time = alpha/(2*pi) log(mu_frag^2/mu_0^2) >= 0.
// Resolves positive log-loss jumps, exponentiates their real-minus-virtual
// generator, and keeps the unresolved/no-resolved-emission endpoint weight.
// LL continuum accuracy requires soft-cutoff/grid convergence. Not NLL.
class LLNonSinglet {
public:
  explicit LLNonSinglet(Real time, const LLGrid& grid = {});
  Real time() const { return time_; }
  const LLGrid& grid() const { return grid_; }
  const LLDiagnostics& diagnostics() const { return diagnostics_; }
  const std::vector<Real>& weights() const { return weights_; }
  // Lattice P0 generator whose exponential gives weights().
  const std::vector<Real>& generator_weights() const { return generator_; }
  Real z(std::size_t index) const;
  Real apply(const Function& measurement) const;
  Real moment(Real N) const;       // numerical, N>=1 (bounded small-z weight)
  Real exact_moment(Real N) const; // continuum exp(time*P0(N)), real N>0

private:
  Real time_;
  LLGrid grid_;
  LLDiagnostics diagnostics_;
  std::vector<Real> weights_,generator_;
};

struct LLBoundaryDiagnostics {
  Real initial_number_residual = 0;
  Real evolved_number_residual = 0;
  Real imaginary_weight_mass = 0;
  Real wrap_coefficient_bound = 0;
  Real omitted_tail_coefficient_bound = 0;
};

// U_LL x (I + a*d1), with the boundary kept LINEAR and allowed to be signed.
// Provides boundary-generated subleading terms, not the missing NLO kernel.
// At two-leg level retain at most one d1 insertion across both muons.
class LLBoundaryNonSinglet {
public:
  explicit LLBoundaryNonSinglet(const QEDParameters& parameters,
                                const LLGrid& grid = {});
  const QEDParameters& parameters() const { return parameters_; }
  const LLNonSinglet& evolution() const { return evolution_; }
  const LLBoundaryDiagnostics& diagnostics() const { return diagnostics_; }
  const std::vector<Real>& initial_boundary_weights() const { return initial_; }
  const std::vector<Real>& evolved_boundary_weights() const { return evolved_; }
  const std::vector<Real>& first_order_weights() const { return order1_; }
  Real apply(const Function& measurement) const;
  Real boundary_moment(Real N) const; // coefficient U_LL x d1, numerical N>=1
  Real exact_boundary_moment(Real N, const IntegrationOptions& options = {}) const;
  Real first_order_apply(const Function& measurement) const;

private:
  QEDParameters parameters_;
  LLNonSinglet evolution_;
  LLBoundaryDiagnostics diagnostics_;
  std::vector<Real> initial_, evolved_, order1_;
};
} // namespace leptonff
