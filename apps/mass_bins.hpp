#pragma once
#include "reference_config.hpp"
#include <sstream>

namespace leptonff {
inline std::vector<Real> read_mass_edges(const std::string& path) {
  std::ifstream input(path);
  if(!input) throw std::runtime_error("Cannot open mass-bin configuration: "+path);
  std::vector<Real> result;
  std::string line;
  bool seen=false;
  while(std::getline(input,line)) {
    line=trim_config(line.substr(0,line.find('#')));
    if(line.empty()) continue;
    const auto equal=line.find('=');
    if(equal==std::string::npos || trim_config(line.substr(0,equal))!="mass_edges" || seen)
      throw std::invalid_argument("Bin configuration must contain exactly one mass_edges entry");
    seen=true;
    const auto values=trim_config(line.substr(equal+1));
    if(values.empty() || values.back()==',') throw std::invalid_argument("Empty mass-bin edge");
    std::istringstream stream(values);
    std::string value;
    while(std::getline(stream,value,',')) {
      value=trim_config(value);
      std::size_t used=0;
      const Real edge=std::stold(value,&used);
      if(used!=value.size() || !std::isfinite(edge))
        throw std::invalid_argument("Invalid mass-bin edge");
      result.push_back(edge);
    }
  }
  if(!seen) throw std::invalid_argument("Missing mass_edges entry");
  return result;
}
} // namespace leptonff
