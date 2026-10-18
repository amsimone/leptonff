#pragma once

#include "leptonff/evolution.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace leptonff {
inline Real parse_ll_real(const char* text) {
  std::size_t consumed = 0;
  const std::string value(text);
  const Real parsed = std::stold(value,&consumed);
  if (consumed != value.size() || !std::isfinite(parsed))
    throw std::invalid_argument("Invalid LL numerical option");
  return parsed;
}

inline LLGrid read_ll_grid(int argc, char** argv) {
  LLGrid grid;
  if (argc >= 3) {
    const Real bins = parse_ll_real(argv[2]);
    if (bins < 16 || bins > (std::size_t{1} << 21) || std::floor(bins) != bins)
      throw std::invalid_argument("Driver bins must be integers in [16,2^21] for refinement");
    grid.bins = static_cast<std::size_t>(bins);
  }
  if (argc >= 4) grid.soft_log_cutoff = parse_ll_real(argv[3]);
  if (argc >= 5) grid.log_span = parse_ll_real(argv[4]);
  grid.validate();
  return grid;
}
} // namespace leptonff
