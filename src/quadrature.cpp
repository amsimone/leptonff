#include "leptonff/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>

namespace leptonff {
namespace {
struct Interval {
  Real lower, upper, value, error;
  bool operator<(const Interval& rhs) const { return error < rhs.error; }
};

Interval rule(const Function& f, Real lower, Real upper) {
  constexpr Real nodes[] = {
      .991455371120812639206854697526329L,
      .949107912342758524526189684047851L,
      .864864423359769072789712788640926L,
      .741531185599394439863864773280788L,
      .586087235467691130294144838258730L,
      .405845151377397166906606412076961L,
      .207784955007898467600689403773245L, 0};
  constexpr Real kronrod[] = {
      .022935322010529224963732008058970L,
      .063092092629978553290700663189204L,
      .104790010322250183839876322541518L,
      .140653259715525918745189590510238L,
      .169004726639267902826583426598550L,
      .190350578064785409913256402421014L,
      .204432940075298892414161999234649L,
      .209482141084727828012999174891714L};
  constexpr Real gauss[] = {
      .129484966168869693270611432679082L,
      .279705391489276667901467771423780L,
      .381830050505118944950369775488975L,
      .417959183673469387755102040816327L};
  const Real mid = lower + (upper - lower) / 2;
  const Real half = (upper - lower) / 2;
  auto sample = [&](Real x) {
    if (!(x > lower && x < upper))
      throw std::runtime_error("Quadrature exhausted floating-point resolution");
    const Real value = f(x);
    if (!std::isfinite(value))
      throw std::runtime_error("Non-finite quadrature sample");
    return value;
  };
  const Real centre = sample(mid);
  Real k = kronrod[7] * centre;
  Real g = gauss[3] * centre;
  Real magnitude = kronrod[7] * std::abs(centre);
  for (int i = 0; i < 7; ++i) {
    const Real left = sample(mid - half * nodes[i]);
    const Real right = sample(mid + half * nodes[i]);
    k += kronrod[i] * (left + right);
    magnitude += kronrod[i] * (std::abs(left) + std::abs(right));
    if (i % 2 == 1) g += gauss[(i - 1) / 2] * (left + right);
  }
  return {lower, upper, half * k,
          std::max(std::abs(half * (k - g)),
                   50 * std::numeric_limits<Real>::epsilon() * half * magnitude)};
}
} // namespace

Integral integrate(const Function& f, Real lower, Real upper,
                   const IntegrationOptions& options,
                   const std::vector<Real>& breakpoints) {
  if (!std::isfinite(lower) || !std::isfinite(upper) || lower > upper)
    throw std::invalid_argument("Integration requires finite ordered bounds");
  if (!std::isfinite(options.absolute_tolerance) ||
      !std::isfinite(options.relative_tolerance) ||
      options.absolute_tolerance <= 0 || options.relative_tolerance < 0 ||
      options.max_intervals == 0)
    throw std::invalid_argument("Invalid integration tolerances or interval limit");
  std::vector<Real> edges{lower, upper};
  for (const Real point : breakpoints) {
    if (!std::isfinite(point))
      throw std::invalid_argument("Non-finite integration breakpoint");
    if (point > lower && point < upper) edges.push_back(point);
  }
  if (lower == upper) return {};
  std::sort(edges.begin(), edges.end());
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
  if (edges.size() - 1 > options.max_intervals)
    throw std::invalid_argument("Too many integration breakpoints");
  std::priority_queue<Interval> intervals;
  Integral result;
  for (std::size_t i = 1; i < edges.size(); ++i) {
    const auto interval = rule(f, edges[i - 1], edges[i]);
    intervals.push(interval);
    result.value += interval.value;
    result.estimated_error += interval.error;
    result.evaluations += 15;
  }
  while (result.estimated_error >
         std::max(options.absolute_tolerance,
                  options.relative_tolerance * std::abs(result.value))) {
    if (intervals.size() >= options.max_intervals)
      throw std::runtime_error("Quadrature did not converge within max_intervals");
    const auto old = intervals.top();
    intervals.pop();
    const Real mid = old.lower + (old.upper - old.lower) / 2;
    const auto left = rule(f, old.lower, mid);
    const auto right = rule(f, mid, old.upper);
    result.value += left.value + right.value - old.value;
    result.estimated_error = std::max(
        Real{0}, result.estimated_error + left.error + right.error - old.error);
    result.evaluations += 30;
    intervals.push(left);
    intervals.push(right);
  }
  return result;
}
} // namespace leptonff
