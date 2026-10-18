#include "leptonff/mass_spectrum.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
// Only these fields are additive measurements. Couplings and compression
// diagnostics are shared metadata and must never be differenced or summed.
constexpr std::array<Real CoupledDecayMatch::*,14> additive{
  &CoupledDecayMatch::born, &CoupledDecayMatch::ff, &CoupledDecayMatch::charge_flipped_pairs,
  &CoupledDecayMatch::muon_projected, &CoupledDecayMatch::muon_geometry,
  &CoupledDecayMatch::photon_projected, &CoupledDecayMatch::photon_geometry,
  &CoupledDecayMatch::lattice_coefficient1, &CoupledDecayMatch::continuum_coefficient1,
  &CoupledDecayMatch::coefficient2, &CoupledDecayMatch::overlap_correction,
  &CoupledDecayMatch::raw, &CoupledDecayMatch::matched, &CoupledDecayMatch::fixed_order};
void metadata(CoupledDecayMatch& out,const CoupledDecayMatch& reference) {
  out.a_initial=reference.a_initial; out.a_hard=reference.a_hard;
  out.compression_moment_residual=reference.compression_moment_residual;
  out.compression_absolute_norm=reference.compression_absolute_norm;
}
}

CoupledMassSpectrum::CoupledMassSpectrum(Real M,const FiducialCuts& cuts,
    const std::vector<Real>& edges,const QEDXSpaceGrid& grid,const CoupledDecayOptions& options)
    : edges_(edges) {
  // Validate binning before the expensive geometry construction.
  if(edges.size()<2 || edges.size()>257)
    throw std::invalid_argument("Mass spectrum needs between 1 and 256 bins");
  for(std::size_t j=0;j<edges.size();++j)
    if(!std::isfinite(edges[j]) || (j && edges[j]<=edges[j-1]))
      throw std::invalid_argument("Mass edges must be finite and strictly increasing");
  if(edges.front()!=cuts.mass_min || edges.back()!=cuts.mass_max)
    throw std::invalid_argument("Mass edges must exactly partition the fiducial window");
  auto histogram_options=options;
  histogram_options.average_mass_cells=true;
  CoupledFiducialDecay base(M,cuts,grid,histogram_options);
  cumulative_.push_back(std::move(base));
  for(std::size_t j=1;j+1<edges.size() && edges[j]<M;++j)
    cumulative_.push_back(cumulative_.front().restricted_mass_window(edges[j],cuts.mass_max));
}

CoupledMassSpectrumResult CoupledMassSpectrum::match(const QEDXSpaceEvolution& evolution,
    const IntegrationOptions& options) const {
  CoupledMassSpectrumResult result;
  std::vector<CoupledDecayMatch> c(edges_.size());
  for(std::size_t j=0;j<cumulative_.size();++j) c[j]=cumulative_[j].match(evolution,options);
  result.total=c.front();
  metadata(result.sum,result.total);
  for(std::size_t j=0;j+1<edges_.size();++j) {
    CoupledMassBin bin;
    bin.lower=edges_[j]; bin.upper=edges_[j+1];
    metadata(bin.integral,result.total);
    for(auto member:additive) {
      bin.integral.*member=c[j].*member-c[j+1].*member;
      result.sum.*member+=bin.integral.*member;
    }
    result.bins.push_back(bin);
  }
  result.matched_closure=result.sum.matched-result.total.matched;
  result.fixed_order_closure=result.sum.fixed_order-result.total.fixed_order;
  result.coefficient1_closure=result.sum.continuum_coefficient1-result.total.continuum_coefficient1;
  result.coefficient2_closure=result.sum.coefficient2-result.total.coefficient2;
  return result;
}
} // namespace leptonff
