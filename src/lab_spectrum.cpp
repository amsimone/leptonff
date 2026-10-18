#include "leptonff/lab_spectrum.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <fstream>
#include <set>
#include <sstream>

namespace leptonff {
Real LabSelection::accepts(bool neutral,const FourVector& a,const FourVector& b) const {
  if(neutral)return z.accepts(a,b)?1:0;
  return a.pt()>w.muon_pt_min && std::abs(a.eta())<w.muon_eta_max &&
    (w.neutrino_pt_min==0 || b.pt()>w.neutrino_pt_min) &&
    (w.transverse_mass_min==0 || transverse_mass(a,b)>w.transverse_mass_min)?1:0;
}
LabSpectrum LabSpectrum::reference(bool neutral,const WLabCuts& w,const FiducialCuts& z) {
  w.validate();z.validate();
  // These edges partition the benchmark fiducial region, including overflows.
  // Reject changed fiducial inputs instead of silently dropping coverage.
  if(w.muon_pt_min!=25 || std::abs(w.muon_eta_max-2.4L)>1e-12L || w.transverse_mass_min!=0 ||
      z.pt_min!=25 || std::abs(z.eta_max-2.4L)>1e-12L || z.mass_min!=60 || z.mass_max!=120)
    throw std::invalid_argument("Reference spectra require pT25, eta2.4, mT0, mass60-120");
  return configured(neutral,reference_binning(neutral),w,z);
}
std::vector<ObservableBinning> LabSpectrum::reference_binning(bool neutral,unsigned subdivision) {
  if(subdivision<1 || subdivision>1024)throw std::invalid_argument("Invalid bin subdivision");
  const Real inf=std::numeric_limits<Real>::infinity();
  std::vector<ObservableBinning> definitions;
  if(neutral) {
    definitions.push_back({"mass",{60,70,80,85,90,95,100,110,120}});
    definitions.push_back({"minpt",{25,30,35,40,50,inf}});
    definitions.push_back({"maxeta",{0,.6L,1.2L,1.8L,2.4L}});
  }else {
    definitions.push_back({"pt",{25,30,35,40,45,60,inf}});
    definitions.push_back({"mt",{0,40,60,70,80,90,120,inf}});
    definitions.push_back({"abseta",{0,.6L,1.2L,1.8L,2.4L}});
  }
  for(auto& d:definitions) {
    std::vector<Real> refined;
    for(std::size_t k=0;k+1<d.edges.size();++k) {
      const unsigned count=std::isfinite(d.edges[k+1])?subdivision:1;
      for(unsigned j=0;j<count;++j)refined.push_back(d.edges[k]+(j==0?0:(d.edges[k+1]-d.edges[k])*j/count));
    }
    refined.push_back(d.edges.back());d.edges=std::move(refined);
  }
  return definitions;
}
std::vector<ObservableBinning> LabSpectrum::read_binning(const std::string& filename) {
  std::ifstream file(filename);
  if(!file)throw std::invalid_argument("Cannot open binning file: "+filename);
  std::vector<ObservableBinning> out;std::string line;
  while(std::getline(file,line)) {
    line=line.substr(0,line.find('#'));std::istringstream row(line);ObservableBinning d;
    if(!(row>>d.observable))continue;
    std::string value;
    while(row>>value) {
      if(value=="inf")d.edges.push_back(std::numeric_limits<Real>::infinity());
      else {std::size_t used=0;const Real x=std::stold(value,&used);
        if(used!=value.size())throw std::invalid_argument("Invalid bin edge: "+value);
        d.edges.push_back(x);}
    }
    out.push_back(std::move(d));
  }
  return out;
}
LabSpectrum LabSpectrum::configured(bool neutral,const std::vector<ObservableBinning>& definitions,
                                    const WLabCuts& w,const FiducialCuts& z) {
  w.validate();z.validate();
  if(definitions.empty())throw std::invalid_argument("No observable bins supplied");
  LabSpectrum out;out.selections.push_back({w,z});
  const Real inf=std::numeric_limits<Real>::infinity();
  const auto append=[&](const std::string& name,const std::vector<Real>& edges,bool increasing,
                        const std::function<void(LabSelection&,Real)>& change) {
    std::vector<int> indices;
    for(std::size_t k=0;k<edges.size();++k) {
      if(k==(increasing?edges.size()-1:0)){indices.push_back(-1);continue;}
      if(k==(increasing?0:edges.size()-1)){indices.push_back(0);continue;}
      auto s=out.selections[0];change(s,edges[k]);
      indices.push_back(static_cast<int>(out.selections.size()));out.selections.push_back(s);
    }
    for(std::size_t k=0;k+1<edges.size();++k)
      out.bins.push_back({name,edges[k],edges[k+1],indices[increasing?k:k+1],indices[increasing?k+1:k]});
  };
  std::set<std::string> names;
  for(const auto& d:definitions) {
    if(!names.insert(d.observable).second || d.edges.size()<2)throw std::invalid_argument("Duplicate/empty observable");
    for(std::size_t i=0;i<d.edges.size();++i)
      if(std::isnan(d.edges[i]) || d.edges[i]<0 || (i && d.edges[i]<=d.edges[i-1]) ||
         (i+1<d.edges.size() && !std::isfinite(d.edges[i])))throw std::invalid_argument("Invalid histogram edges");
    Real lower=0,upper=0;bool increasing=true;
    std::function<void(LabSelection&,Real)> change;
    if(neutral && d.observable=="mass") {lower=z.mass_min;upper=z.mass_max;change=[](LabSelection& s,Real x){s.z.mass_min=x;};}
    else if(neutral && d.observable=="minpt") {lower=z.pt_min;upper=inf;change=[](LabSelection& s,Real x){s.z.pt_min=x;};}
    else if(neutral && d.observable=="maxeta") {upper=z.eta_max;increasing=false;change=[](LabSelection& s,Real x){s.z.eta_max=x;};}
    else if(!neutral && d.observable=="pt") {lower=w.muon_pt_min;upper=inf;change=[](LabSelection& s,Real x){s.w.muon_pt_min=x;};}
    else if(!neutral && d.observable=="mt") {lower=w.transverse_mass_min;upper=inf;change=[](LabSelection& s,Real x){s.w.transverse_mass_min=x;};}
    else if(!neutral && d.observable=="abseta") {upper=w.muon_eta_max;increasing=false;change=[](LabSelection& s,Real x){s.w.muon_eta_max=x;};}
    else throw std::invalid_argument("Unsupported observable for this process: "+d.observable);
    const auto equivalent=[](Real a,Real b){return a==b || (std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=1e-12L*std::max(Real{1},std::max(std::abs(a),std::abs(b))));};
    if(!equivalent(d.edges.front(),lower) || !equivalent(d.edges.back(),upper))throw std::invalid_argument("Bins must cover the full fiducial support: "+d.observable);
    append(d.observable,d.edges,increasing,change);
  }
  return out;
}
Real LabSpectrum::bin(std::size_t k,const std::vector<Real>& values) const {
  if(values.size()!=selections.size())throw std::invalid_argument("Incomplete cumulative spectrum");
  const auto& b=bins.at(k);
  return (b.positive<0?0:values.at(b.positive))-(b.negative<0?0:values.at(b.negative));
}

MeasurementPlan::MeasurementPlan(bool neutral,LabSpectrum spectrum):neutral_(neutral),spectrum_(std::move(spectrum)) {
  if(spectrum_.selections.empty())throw std::invalid_argument("Empty measurement plan");
  for(const auto& s:spectrum_.selections){s.w.validate();s.z.validate();}
  const auto size=spectrum_.selections.size();representatives_.resize(size);
  for(const auto& b:spectrum_.bins)if(b.positive < -1 || b.negative < -1 ||
      b.positive>=static_cast<int>(size) || b.negative>=static_cast<int>(size))throw std::invalid_argument("Invalid measurement bin index");
  for(std::size_t i=0;i<size;++i) {
    const auto& a=spectrum_.selections[i];std::size_t best=i;
    for(std::size_t j=0;j<size;++j) {
      const auto& b=spectrum_.selections[j];
      const bool equal=neutral?(a.z.pt_min==b.z.pt_min && a.z.mass_min==b.z.mass_min && a.z.mass_max==b.z.mass_max):
        (a.w.muon_pt_min==b.w.muon_pt_min && a.w.neutrino_pt_min==b.w.neutrino_pt_min && a.w.transverse_mass_min==b.w.transverse_mass_min);
      if(equal && (neutral?b.z.eta_max>spectrum_.selections[best].z.eta_max:b.w.muon_eta_max>spectrum_.selections[best].w.muon_eta_max))best=j;
    }
    representatives_[i]=best;
  }
}
bool MeasurementPlan::accepts_eta(std::size_t i,const FourVector& a,const FourVector& b) const {
  const auto& s=spectrum_.selections.at(i);
  return neutral_?(a.pt()>0 && b.pt()>0 && std::abs(a.eta())<s.z.eta_max && std::abs(b.eta())<s.z.eta_max):
    (a.pt()>0 && std::abs(a.eta())<s.w.muon_eta_max);
}
} // namespace leptonff
