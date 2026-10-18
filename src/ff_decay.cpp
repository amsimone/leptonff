#include "leptonff/ff_decay.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff {
XSpaceFFDecayMatch decay_xspace_ff_match(Real mass,const FiducialCuts& cuts,
    const QEDXSpaceEvolution& evolution,const IntegrationOptions& options) {
  cuts.validate();
  if (!std::isfinite(mass) || mass<=0) throw std::invalid_argument("Invalid parent mass");
  const Real minimum=std::max(2*cuts.pt_min/mass,std::pow(cuts.mass_min/mass,2));
  if (minimum<=0) throw std::invalid_argument("Inclusive pair-density reference requires a positive pT or mass lower cut");
  const auto& grid=evolution.grid();
  if(minimum<evolution.z(grid.bins-1))
    throw std::invalid_argument("Fiducial support extends below the coupled grid");
  auto undo_tilt = [&](const std::vector<Real>& weights) {
    std::vector<Real> result(grid.bins,0);
    for(std::size_t k=0;k<grid.bins;++k) {
      const Real z=evolution.z(k);
      if(z<minimum) break;
      result[k]=weights[k]/z;
    }
    return result;
  };
  const auto same=undo_tilt(evolution.weights()[0][0]);
  const auto opposite=undo_tilt(evolution.weights()[0][1]);
  const auto first=undo_tilt(evolution.coefficient1()[0][0]);
  std::vector<Real> identity(grid.bins,0); identity[0]=1;
  XSpaceFFDecayMatch result;
  result.born=isotropic_acceptance(mass,1,1,cuts);
  result.charge_flipped_pairs=decay_lattice_acceptance(mass,cuts,grid.spacing(),opposite,opposite);
  result.resummed=decay_lattice_acceptance(mass,cuts,grid.spacing(),same,same)+result.charge_flipped_pairs;
  result.lattice_coefficient1=2*decay_lattice_acceptance(mass,cuts,grid.spacing(),first,identity);
  const auto& p=evolution.parameters();
  const QEDParameters q{p.alpha_msbar_at_mu0,p.muon_mass,p.mu_0,p.mu_frag};
  result.continuum_coefficient1=decay_ff_coefficient_order1(mass,cuts,q,options);
  result.additive_correction=q.expansion_parameter()*(result.continuum_coefficient1-result.lattice_coefficient1);
  result.matched=result.resummed+result.additive_correction;
  return result;
}
} // namespace leptonff
