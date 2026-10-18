#pragma once
#include "leptonff/vector_decay.hpp"
#include <cmath>

namespace leptonff::detail {
// Integral from zero to exp(-u). Separate soft/tail branches avoid
// cancellation; the whole-plus endpoint is handled by the caller.
inline Real hard_primitive(Real u,bool w_extra) {
  const Real z=std::exp(-u),t=-std::expm1(-u);
  const Real l=u<.5L?std::log(t):std::log1p(-z);
  if(w_extra)return 2*l+4*z/3+z*z/3;
  if(u>=1)return (2*l+z+z*z/2)*u-2*real_dilogarithm(z)+4*z+3*z*z/4
                 -2*l*l+(3-2*z-z*z)*l;
  const Real lz=2*real_dilogarithm(t)+(z+z*z/2)*u+z+z*z/4;
  const Real ls=-2*l*l+(4*t-t*t)*l-4*t+t*t/2;
  return lz+ls-std::pow(std::acos(-1.L),2)/3+3.5L;
}
inline Real projected_hard_primitive(Real z,bool charged,Real scale_log) {
  if(z==0)return 0;
  const Real u=-std::log(z);
  return hard_primitive(u,false)+(charged?hard_primitive(u,true):0)
    +scale_log*(-2*std::log1p(-z)-z-z*z/2);
}
} // namespace leptonff::detail
