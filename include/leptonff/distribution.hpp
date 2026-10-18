#pragma once

#include "leptonff/quadrature.hpp"

namespace leptonff {
// delta * delta(1-z) + regular(z) + sum_i [plus_i(z)]_+.
// The plus prescription acts on the WHOLE kernel, including its numerator.
class Distribution {
public:
  Distribution(Real delta = 0, Function regular = {},
               std::vector<Function> plus = {}, Real moment_min_exclusive = 0);
  Real apply(const Function& measurement, Real lower = 0,
             const IntegrationOptions& options = {},
             const std::vector<Real>& breakpoints = {}) const;
  Real moment(Real N, const IntegrationOptions& options = {}) const;
  // Integral over [lower,upper); the final bin with upper=1 includes the
  // endpoint distribution. Known edges are supplied to the quadrature.
  Real integral(Real lower, Real upper,
                const IntegrationOptions& options = {}) const;
  // (D tensor f)(x) = integral_x^1 dz/z D(z) f(x/z), including the
  // endpoint and whole-plus subtraction. f must be smooth on [x,1].
  Real convolve(const Function& f, Real x,
                const IntegrationOptions& options = {}) const;
  // Interior real-emission density regular(z)+sum plus_i(z), 0<z<1.
  // Does not include delta/endpoint subtraction; not a complete FF by itself.
  Real resolved_density(Real z) const;
  Distribution scaled(Real factor) const;
  friend Distribution operator+(const Distribution&, const Distribution&);

private:
  Real delta_;
  Function regular_;
  std::vector<Function> plus_;
  Real moment_min_exclusive_;
};

using PairMeasurement = std::function<Real(Real, Real)>;
using InnerBreakpoints = std::function<std::vector<Real>(Real)>;

// Full tensor product, including all delta/plus cross and double-endpoint terms.
// For discontinuous cuts supply outer and conditional inner breakpoints.
Real apply_pair(const Distribution& first, const Distribution& second,
                const PairMeasurement& measurement,
                const IntegrationOptions& options = {},
                const std::vector<Real>& outer_breakpoints = {},
                const InnerBreakpoints& inner_breakpoints = {});

// I x I + a (D1_first x I + I x D1_second), with no spurious a^2 term.
Real apply_pair_order1(const Distribution& first_coefficient,
                       const Distribution& second_coefficient, Real a,
                       const PairMeasurement& measurement,
                       const IntegrationOptions& options = {},
                       const std::vector<Real>& first_breakpoints = {},
                       const std::vector<Real>& second_breakpoints = {});
} // namespace leptonff
