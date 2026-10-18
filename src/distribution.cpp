#include "leptonff/distribution.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace leptonff {
Distribution::Distribution(Real delta, Function regular,
                           std::vector<Function> plus, Real moment_min_exclusive)
    : delta_(delta), regular_(std::move(regular)), plus_(std::move(plus)),
      moment_min_exclusive_(moment_min_exclusive) {
  if (!std::isfinite(delta_) || !std::isfinite(moment_min_exclusive_))
    throw std::invalid_argument("Non-finite distribution parameter");
  for (const auto& kernel : plus_)
    if (!kernel) throw std::invalid_argument("Empty plus kernel");
}

Real Distribution::apply(const Function& measurement, Real lower,
                          const IntegrationOptions& options,
                          const std::vector<Real>& breakpoints) const {
  // lower=1 is excluded: a plus distribution has no finite action on a
  // shrinking interval containing only its singular endpoint.
  if (!std::isfinite(lower) || lower < 0 || lower >= 1)
    throw std::invalid_argument("Distribution lower bound must be in [0,1)");
  const Real endpoint = measurement(1);
  if (!std::isfinite(endpoint))
    throw std::invalid_argument("Measurement must be finite at z=1");
  Real value = delta_ * endpoint;
  if (regular_)
    value += integrate([&](Real z) { return regular_(z) * measurement(z); },
                       lower, 1, options, breakpoints).value;
  for (const auto& kernel : plus_) {
    value += integrate([&](Real z) {
      return kernel(z) * (measurement(z) - endpoint);
    }, lower, 1, options, breakpoints).value;
    if (lower > 0 && endpoint != 0)
      value -= endpoint * integrate(kernel, 0, lower, options).value;
  }
  if (!std::isfinite(value))
    throw std::runtime_error("Non-finite distribution action");
  return value;
}

Real Distribution::moment(Real N, const IntegrationOptions& options) const {
  if (!std::isfinite(N) || N <= moment_min_exclusive_)
    throw std::invalid_argument("Moment is outside the distribution's domain");
  return apply([N](Real z) { return std::pow(z, N - 1); }, 0, options);
}

Real Distribution::integral(Real lower, Real upper,
                            const IntegrationOptions& options) const {
  if (!std::isfinite(lower) || !std::isfinite(upper) ||
      lower < 0 || lower >= upper || upper > 1)
    throw std::invalid_argument("Distribution bin requires 0<=lower<upper<=1");
  return apply([upper](Real z) { return z <= upper ? 1.0L : 0.0L; },
               lower, options, upper < 1 ? std::vector<Real>{upper} : std::vector<Real>{});
}

Real Distribution::convolve(const Function& f, Real x,
                            const IntegrationOptions& options) const {
  if (!f || !std::isfinite(x) || x <= 0 || x >= 1)
    throw std::invalid_argument("Convolution requires a function and 0<x<1");
  return apply([&](Real z) { return f(x/z)/z; }, x, options);
}

Distribution Distribution::scaled(Real factor) const {
  if (!std::isfinite(factor)) throw std::invalid_argument("Non-finite scaling");
  if (factor == 0) return {};
  Function regular;
  if (regular_) regular = [f = regular_, factor](Real z) { return factor * f(z); };
  std::vector<Function> plus;
  for (const auto& f : plus_)
    plus.emplace_back([f, factor](Real z) { return factor * f(z); });
  return {factor * delta_, regular, plus, moment_min_exclusive_};
}

Real Distribution::resolved_density(Real z) const {
  if (!std::isfinite(z) || z <= 0 || z >= 1)
    throw std::invalid_argument("Resolved density requires 0<z<1");
  Real result = regular_ ? regular_(z) : 0;
  for (const auto& kernel : plus_) result += kernel(z);
  if (!std::isfinite(result)) throw std::runtime_error("Non-finite resolved density");
  return result;
}

Distribution operator+(const Distribution& lhs, const Distribution& rhs) {
  Function regular;
  if (lhs.regular_ || rhs.regular_)
    regular = [l = lhs.regular_, r = rhs.regular_](Real z) {
      return (l ? l(z) : 0) + (r ? r(z) : 0);
    };
  auto plus = lhs.plus_;
  plus.insert(plus.end(), rhs.plus_.begin(), rhs.plus_.end());
  return {lhs.delta_ + rhs.delta_, regular, plus,
          std::max(lhs.moment_min_exclusive_, rhs.moment_min_exclusive_)};
}

Real apply_pair(const Distribution& first, const Distribution& second,
                const PairMeasurement& measurement,
                const IntegrationOptions& options,
                const std::vector<Real>& outer_breakpoints,
                const InnerBreakpoints& inner_breakpoints) {
  return first.apply([&](Real z1) {
    return second.apply([&](Real z2) { return measurement(z1, z2); }, 0, options,
                        inner_breakpoints ? inner_breakpoints(z1) : std::vector<Real>{});
  }, 0, options, outer_breakpoints);
}

Real apply_pair_order1(const Distribution& first_coefficient,
                       const Distribution& second_coefficient, Real a,
                       const PairMeasurement& measurement,
                       const IntegrationOptions& options,
                       const std::vector<Real>& first_breakpoints,
                       const std::vector<Real>& second_breakpoints) {
  if (!std::isfinite(a)) throw std::invalid_argument("Non-finite coupling");
  return measurement(1, 1) + a * (
      first_coefficient.apply([&](Real z) { return measurement(z, 1); },
                              0, options, first_breakpoints) +
      second_coefficient.apply([&](Real z) { return measurement(1, z); },
                               0, options, second_breakpoints));
}
} // namespace leptonff
