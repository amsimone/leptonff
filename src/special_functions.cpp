#include "leptonff/special_functions.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

namespace leptonff {
Real real_dilogarithm(Real x) {
  if (!std::isfinite(x) || x > 1)
    throw std::invalid_argument("Real dilogarithm requires finite x<=1");
  const Real zeta2 = std::pow(std::acos(-1.0L),2)/6;
  if (x == 0) return 0;
  if (x == 1) return zeta2;
  if (x < 0)
    return -real_dilogarithm(x/(x-1))-std::pow(std::log1p(-x),2)/2;
  if (x > .5L)
    return zeta2-std::log(x)*std::log1p(-x)-real_dilogarithm(1-x);
  Real sum = 0, power = x;
  for (unsigned k = 1; k < 1024; ++k) {
    const Real term = power/(Real(k)*k);
    sum += term;
    power *= x;
    if (std::abs(term) < std::numeric_limits<Real>::epsilon()*std::abs(sum)/2)
      return sum;
  }
  throw std::runtime_error("Dilogarithm series did not converge");
}

} // namespace leptonff
