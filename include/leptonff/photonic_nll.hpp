#pragma once
#include "leptonff/evolution.hpp"

namespace leptonff {
enum class PhotonicAccuracy { ll, nll };

// Published time-like photonic diagram class G has G(N=1) != 0. Keep that
// source intact; this operator describes the explicitly conditioned shape
// E(N)/E(1), consistent with the unitary decay/branching-ratio convention.
Real photonic_source_number_moment();
Distribution unitary_photonic_nlo_kernel();

struct PhotonicNLLDiagnostics {
  Real source_number_moment = 0, raw_inclusive_factor = 1;
  Real kernel_number_residual = 0, evolved_number_residual = 0;
  Real imaginary_weight_mass = 0, wrap_bound = 0, omitted_tail_bound = 0;
};

// U_LL [1 + a*d1 + a^2*L*(G-G(1)*delta)]. Exactly one subleading
// insertion: excludes d1*G, G*G and two-leg products of finite insertions.
// Fixed alpha, no real pairs or fermion-loop running. Signed weights retained.
class PhotonicNLL {
public:
  PhotonicNLL(const QEDParameters&,const LLGrid& = {});
  explicit PhotonicNLL(const LLBoundaryNonSinglet&);
  const LLBoundaryNonSinglet& baseline() const { return baseline_; }
  const PhotonicNLLDiagnostics& diagnostics() const { return diagnostics_; }
  const std::vector<Real>& kernel_weights() const { return kernel_; }
  const std::vector<Real>& evolved_kernel_weights() const { return evolved_; }
  // Coefficient multiplying a: U*d1 + a*L*U*(G-G(1)*delta).
  const std::vector<Real>& finite_insertion_weights() const { return insertion_; }
  const std::vector<Real>& coefficient2_weights() const { return order2_; }
  Real apply(const Function&) const;
  Real moment(Real N) const;
  Real exact_moment(Real N,const IntegrationOptions& = {}) const;
private:
  LLBoundaryNonSinglet baseline_;
  PhotonicNLLDiagnostics diagnostics_;
  std::vector<Real> kernel_,evolved_,insertion_,order2_;
};
} // namespace leptonff
