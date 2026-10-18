#pragma once

#include "leptonff/quadrature.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

inline void close(leptonff::Real actual, leptonff::Real expected,
                  leptonff::Real tolerance, const std::string& name) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
    std::cerr.precision(18);
    std::cerr << name << ": actual=" << actual << " expected=" << expected << '\n';
    throw std::runtime_error(name);
  }
}

template <typename Exception = std::exception, typename F>
void throws(F&& action, const std::string& name) {
  try { action(); } catch (const Exception&) { return; }
  throw std::runtime_error(name + ": expected exception");
}
