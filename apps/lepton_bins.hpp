#pragma once
#include "reference_config.hpp"
#include <sstream>

namespace leptonff {
struct LeptonBins {std::vector<Real> pt,eta;};
inline LeptonBins read_lepton_bins(const std::string& path) {
  std::ifstream input(path);
  if(!input)throw std::runtime_error("Cannot open lepton-bin config: "+path);
  LeptonBins result;
  std::set<std::string> seen;
  std::string line;
  while(std::getline(input,line)) {
    line=trim_config(line.substr(0,line.find('#')));if(line.empty())continue;
    const auto equal=line.find('=');
    if(equal==std::string::npos)throw std::invalid_argument("Expected key = values in lepton-bin config");
    const auto key=trim_config(line.substr(0,equal)),values=trim_config(line.substr(equal+1));
    if((key!="pt_edges" && key!="eta_edges") || !seen.insert(key).second || values.empty() || values.back()==',')
      throw std::invalid_argument("Expected unique pt_edges and eta_edges entries");
    auto& edges=key=="pt_edges"?result.pt:result.eta;
    std::istringstream stream(values);std::string value;
    while(std::getline(stream,value,',')) {
      value=trim_config(value);std::size_t used=0;const Real x=std::stold(value,&used);
      if(used!=value.size() || !std::isfinite(x))throw std::invalid_argument("Invalid lepton-bin edge");
      edges.push_back(x);
    }
  }
  if(seen.size()!=2)throw std::invalid_argument("Both pt_edges and eta_edges are required");
  return result;
}
}
