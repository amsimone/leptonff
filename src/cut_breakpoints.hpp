#pragma once
#include "leptonff/quadrature.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace leptonff::detail {
// Algebraically equivalent angular contacts can round to neighbouring values
// (in particular after folding c and -c). Do not create intervals too narrow
// to contain all 15 quadrature nodes. This acts only on generated finite-cut
// contacts, never on a physical integration endpoint or a source distribution.
inline std::vector<Real> resolved_cut_breakpoints(std::vector<Real> points,Real lower,Real upper) {
  std::sort(points.begin(),points.end());
  std::vector<Real> result;
  const Real margin=1024*std::numeric_limits<Real>::epsilon()
      *std::max({Real{1},std::abs(lower),std::abs(upper)});
  Real previous=lower;
  for(Real point:points) if(point-previous>margin && upper-point>margin) {
    result.push_back(point);previous=point;
  }
  return result;
}
} // namespace leptonff::detail
