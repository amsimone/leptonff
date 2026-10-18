#pragma once

#include "leptonff/decay.hpp"

#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace leptonff {
struct ReferenceInput {
  Real mass = 91.1876L;
  QEDParameters qed;
  FiducialCuts cuts;
  IntegrationOptions integration;
};

inline std::string trim_config(const std::string& text) {
  const auto first = text.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) return {};
  return text.substr(first, text.find_last_not_of(" \t\r\n")-first+1);
}

inline void read_reference_config(const std::string& path, ReferenceInput& result) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open configuration: " + path);
  const std::map<std::string, Real*> fields{
      {"parent_mass", &result.mass}, {"muon_mass", &result.qed.muon_mass},
      {"mu_0", &result.qed.mu_0}, {"mu_frag", &result.qed.mu_frag},
      {"alpha", &result.qed.alpha}, {"pt_min", &result.cuts.pt_min},
      {"eta_max", &result.cuts.eta_max}, {"mass_min", &result.cuts.mass_min},
      {"mass_max", &result.cuts.mass_max},
      {"absolute_tolerance", &result.integration.absolute_tolerance},
      {"relative_tolerance", &result.integration.relative_tolerance}};
  std::set<std::string> seen;
  std::string line;
  std::size_t number = 0;
  while (std::getline(input, line)) {
    ++number;
    line = trim_config(line.substr(0, line.find('#')));
    if (line.empty()) continue;
    const auto equal = line.find('=');
    if (equal == std::string::npos)
      throw std::runtime_error("Expected key = value at line " + std::to_string(number));
    const auto key = trim_config(line.substr(0, equal));
    const auto value = trim_config(line.substr(equal+1));
    const auto field = fields.find(key);
    if (field == fields.end()) throw std::runtime_error("Unknown configuration key: " + key);
    if (!seen.insert(key).second) throw std::runtime_error("Duplicate configuration key: " + key);
    std::size_t consumed = 0;
    const Real parsed = std::stold(value, &consumed);
    if (consumed != value.size() || !std::isfinite(parsed))
      throw std::runtime_error("Invalid numeric value for " + key);
    *field->second = parsed;
  }
  if (seen.size() != fields.size())
    throw std::runtime_error("Configuration must explicitly set every reference parameter");
}
} // namespace leptonff
