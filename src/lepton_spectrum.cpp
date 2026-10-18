#include "leptonff/lepton_spectrum.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
constexpr std::array<Real CoupledDecayMatch::*,14> additive{
  &CoupledDecayMatch::born,&CoupledDecayMatch::ff,&CoupledDecayMatch::charge_flipped_pairs,
  &CoupledDecayMatch::muon_projected,&CoupledDecayMatch::muon_geometry,
  &CoupledDecayMatch::photon_projected,&CoupledDecayMatch::photon_geometry,
  &CoupledDecayMatch::lattice_coefficient1,&CoupledDecayMatch::continuum_coefficient1,
  &CoupledDecayMatch::coefficient2,&CoupledDecayMatch::overlap_correction,
  &CoupledDecayMatch::raw,&CoupledDecayMatch::matched,&CoupledDecayMatch::fixed_order};
void metadata(CoupledDecayMatch& out,const CoupledDecayMatch& reference) {
  out.a_initial=reference.a_initial;out.a_hard=reference.a_hard;
  out.compression_moment_residual=reference.compression_moment_residual;
  out.compression_absolute_norm=reference.compression_absolute_norm;
}
}

CoupledDecayOptions lepton_spectrum_options(Real mass,const FiducialCuts& cuts,
    const std::vector<Real>& pt_edges,const std::vector<Real>& eta_edges,const CoupledDecayOptions& base) {
  cuts.validate();base.validate();
  if(!std::isfinite(mass) || mass<=0) throw std::invalid_argument("Invalid spectrum mass");
  auto result=base;result.exact_endpoint_strips=true;result.average_mass_cells=true;
  const auto contact=[&](Real f) {if(f>0 && f<1) result.ff_fraction_breaks.push_back(f);};
  for(Real t:pt_edges) {
    if(!std::isfinite(t)) throw std::invalid_argument("Non-finite pT edge");
    contact(2*t/mass);
  }
  for(Real t:eta_edges) {
    if(!std::isfinite(t)) throw std::invalid_argument("Non-finite eta edge");
    if(std::abs(t)<=cuts.eta_max && std::tanh(std::abs(t))<1)
      contact(2*cuts.pt_min*std::cosh(t)/mass);
  }
  std::sort(result.ff_fraction_breaks.begin(),result.ff_fraction_breaks.end());
  result.ff_fraction_breaks.erase(std::unique(result.ff_fraction_breaks.begin(),result.ff_fraction_breaks.end()),
      result.ff_fraction_breaks.end());
  result.validate();return result;
}

CoupledLeptonSpectrum::CoupledLeptonSpectrum(Real M,const FiducialCuts& cuts,LeptonObservable observable,
    const std::vector<Real>& edges,const QEDXSpaceGrid& grid,const CoupledDecayOptions& options)
    :observable_(observable),edges_(edges) {
  if(edges.size()<2 || edges.size()>257)
    throw std::invalid_argument("Single-lepton spectrum needs 1 to 256 bins");
  cuts.validate();
  for(std::size_t j=0;j<edges.size();++j)
    if(!std::isfinite(edges[j]) || (j && edges[j]<=edges[j-1]))
      throw std::invalid_argument("Single-lepton edges must be finite and strictly increasing");
  if(observable==LeptonObservable::transverse_momentum) {
    if(edges.front()!=cuts.pt_min || edges.back()<M/2)
      throw std::invalid_argument("pT edges must start at pt_min and cover the kinematic endpoint M/2");
    thresholds_=edges;
  } else if(observable==LeptonObservable::pseudorapidity) {
    if(edges.front()!=-cuts.eta_max || edges.back()!=cuts.eta_max)
      throw std::invalid_argument("Eta edges must partition [-eta_max,+eta_max]");
    thresholds_.push_back(0);
    for(Real edge:edges) thresholds_.push_back(std::abs(edge));
    std::sort(thresholds_.begin(),thresholds_.end());
    thresholds_.erase(std::unique(thresholds_.begin(),thresholds_.end()),thresholds_.end());
  } else throw std::invalid_argument("Unknown single-lepton observable");
  auto histogram_options=options;histogram_options.average_mass_cells=true;histogram_options.exact_endpoint_strips=true;
  cumulative_.emplace_back(M,cuts,grid,histogram_options);
  for(Real t:thresholds_) {
    if(observable==LeptonObservable::transverse_momentum) {
      if(t<=cuts.pt_min) {indices_.push_back(-2);continue;}
      if(t>=M/2) {indices_.push_back(-1);continue;}
    } else {
      if(t==0) {indices_.push_back(-1);continue;}
      if(t>=cuts.eta_max) {indices_.push_back(-2);continue;}
    }
    const auto cut=observable==LeptonObservable::transverse_momentum
        ? SingleLeptonCut::pt_above : SingleLeptonCut::absolute_eta_below;
    cumulative_.push_back(cumulative_.front().with_lepton_selection({cut,t}));
    indices_.push_back(static_cast<int>(cumulative_.size()-1));
  }
}

CoupledMassSpectrumResult CoupledLeptonSpectrum::match(const QEDXSpaceEvolution& e,
    const IntegrationOptions& options) const {
  CoupledMassSpectrumResult result;
  std::vector<CoupledDecayMatch> c;
  for(const auto& op:cumulative_) c.push_back(op.match(e,options));
  result.total=c.front();metadata(result.sum,result.total);
  const CoupledDecayMatch zero;
  const auto value=[&](Real t)->const CoupledDecayMatch& {
    const auto index=static_cast<std::size_t>(std::lower_bound(thresholds_.begin(),thresholds_.end(),t)-thresholds_.begin());
    const int location=indices_[index];
    return location==-1 ? zero : (location==-2 ? result.total : c[static_cast<std::size_t>(location)]);
  };
  for(std::size_t j=0;j+1<edges_.size();++j) {
    CoupledMassBin bin;bin.lower=edges_[j];bin.upper=edges_[j+1];metadata(bin.integral,result.total);
    for(auto member:additive) {
      Real x=0;
      if(observable_==LeptonObservable::transverse_momentum)
        x=2*(value(bin.lower).*member-value(bin.upper).*member);
      else if(bin.upper<=0) x=value(-bin.lower).*member-value(-bin.upper).*member;
      else if(bin.lower>=0) x=value(bin.upper).*member-value(bin.lower).*member;
      else x=value(-bin.lower).*member+value(bin.upper).*member;
      bin.integral.*member=x;result.sum.*member+=x;
    }
    result.bins.push_back(bin);
  }
  result.matched_closure=result.sum.matched-2*result.total.matched;
  result.fixed_order_closure=result.sum.fixed_order-2*result.total.fixed_order;
  result.coefficient1_closure=result.sum.continuum_coefficient1-2*result.total.continuum_coefficient1;
  result.coefficient2_closure=result.sum.coefficient2-2*result.total.coefficient2;
  return result;
}
} // namespace leptonff
