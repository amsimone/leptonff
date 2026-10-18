#pragma once

#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace leptonff {
using Real = long double;
#ifdef LEPTONFF_REAL_MANTISSA
static_assert(std::numeric_limits<Real>::digits==LEPTONFF_REAL_MANTISSA,
  "leptonff native C++ long-double ABI mismatch: use a matching compiler or the double-based C interface");
#endif
using Function = std::function<Real(Real)>;

struct IntegrationOptions {
  Real absolute_tolerance = 1e-10L;
  Real relative_tolerance = 1e-10L;
  std::size_t max_intervals = 8192;
};

struct Integral {
  Real value = 0;
  Real estimated_error = 0;
  std::size_t evaluations = 0;
};

// Adaptive embedded Gauss(7)-Kronrod(15). Breakpoints must be supplied for
// known discontinuities; no adaptive rule can guarantee finding a narrow bin.
// Non-finite samples and failure to reach tolerance throw exceptions.
Integral integrate(const Function& f, Real lower, Real upper,
                   const IntegrationOptions& options = {},
                   const std::vector<Real>& breakpoints = {});
} // namespace leptonff
