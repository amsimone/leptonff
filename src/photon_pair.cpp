#include "leptonff/photon_pair.hpp"
#include "cut_breakpoints.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
void positive(Real x) {
  if (!std::isfinite(x) || x <= 0)
    throw std::invalid_argument("Photon-pair scales must be positive and finite");
}
void fraction(Real f) {
  if (!std::isfinite(f) || f < 0 || f > 1)
    throw std::invalid_argument("Photon-pair fractions must be in [0,1]");
}
Real bound(Real energy, const FiducialCuts& cuts) {
  if (energy <= cuts.pt_min) return 0;
  const Real r = cuts.pt_min/energy;
  return std::min(std::tanh(cuts.eta_max),std::sqrt((1-r)*(1+r)));
}
Real acceptance(Real M, Real y, Real t, Real fg, Real fm,
                const FiducialCuts& cuts,const SingleLeptonSelection& selection = {}) {
  const Real m = M*std::sqrt(y*fg*fm);
  if (m <= cuts.mass_min || m >= cuts.mass_max) return 0;
  const Real xg = y+t, xm = 1-t;
  const Real bg = bound(M*xg*fg/2,cuts), bm = bound(M*xm*fm/2,cuts);
  if (bg == 0 || bm == 0) return 0;
  if(selection.cut!=SingleLeptonCut::none)
    return rotational_lepton_acceptance(M*xg*fg/2,M*xm*fm/2,
        std::clamp(1-2*y/(xg*xm),Real{-1},Real{1}),cuts,selection);
  return rotational_band_acceptance(bg,bm,
      std::clamp(1-2*y/(xg*xm),Real{-1},Real{1}));
}
} // namespace

Distribution vector_photon_pair_mass_hard_coefficient(Real M, Real mu) {
  positive(M); positive(mu);
  const Real L = 2*(std::log(M)-std::log(mu));
  return {0,[L](Real y) {
    const Real v = 1-y, g = (1+v*v)/y;
    return g*(L+2*std::log1p(-y)+std::log(y)+1)
        -4*v/y+v*v/(2*y);
  },{},1};
}

Real massive_vector_photon_pair_spectrum(Real M, Real m, Real y) {
  positive(M); positive(m);
  if (m >= M/2 || !std::isfinite(y) || y <= 0 || y >= 1)
    throw std::invalid_argument("Invalid massive photon-pair spectrum point");
  if (y >= 1-2*m/M) return 0;
  const Real r = std::pow(m/M,2), v = 1-y;
  if (r == 0) throw std::invalid_argument("Mass ratio underflows");
  const Real root = std::sqrt((v-2*m/M)*(v+2*m/M));
  const Real tmax = y*(v-2*r+root)/(2*(y+r));
  const Real tmin = 2*r*y/(v-2*r+root);
  const Real pole = (1+v*v)/y-4*r-8*r*r/y;
  const Real constant = -(2+4*r)/y-2*r*(1+2*r)/(y*y);
  return pole*std::log(tmax/tmin)+constant*(tmax-tmin)
      +(tmax-tmin)*(tmax+tmin)/(2*y)
      +2*r*(1+2*r)*(1/tmax-1/tmin);
}

Real vector_photon_pair_acceptance(Real M, Real y, Real t,
    Real fg, Real fm, const FiducialCuts& cuts,const SingleLeptonSelection& selection) {
  positive(M); cuts.validate(); fraction(fg); fraction(fm);
  selection.validate();
  if (!std::isfinite(y) || y <= 0 || y >= 1 || !std::isfinite(t) ||
      t < 0 || t > 1-y)
    throw std::invalid_argument("Invalid photon-pair invariants");
  return acceptance(M,y,t,fg,fm,cuts,selection);
}

PhotonPairHardCoefficient fiducial_vector_photon_pair_hard(Real M,
    Real mu, const FiducialCuts& cuts, Real fg, Real fm,
    const HardDecayOptions& options,const SingleLeptonSelection& selection) {
  positive(M); positive(mu); cuts.validate(); fraction(fg); fraction(fm);
  options.validate();
  selection.validate();
  if (cuts.mass_min <= 0)
    throw std::invalid_argument("Photon-pair hard action needs a positive pair-mass cut");
  PhotonPairHardCoefficient result;
  if (fg == 0 || fm == 0) return result;
  const Real lower = std::pow(cuts.mass_min/M,2)/(fg*fm);
  const Real upper = std::min(Real{1},std::pow(cuts.mass_max/M,2)/(fg*fm));
  if (lower >= upper || 2*cuts.pt_min/M >= std::min(fg,fm)) return result;
  const auto B = [&](Real y) { return isotropic_acceptance(M,y*fg,fm,cuts,selection); };
  std::vector<Real> breaks{2*cuts.pt_min/(M*fg),fm/fg,
      2*cuts.pt_min/M*(1/fg+1/fm)-1};
  if (cuts.pt_min > 0 && std::tanh(cuts.eta_max) < 1)
    breaks.push_back(2*cuts.pt_min*std::cosh(cuts.eta_max)/(M*fg));
  std::vector<Real> energies{cuts.pt_min};
  if(selection.cut==SingleLeptonCut::pt_above) {
    energies.push_back(selection.threshold);
    if(std::tanh(cuts.eta_max)<1) energies.push_back(selection.threshold*std::cosh(cuts.eta_max));
    if(selection.threshold>0) breaks.push_back(fm/fg*cuts.pt_min/selection.threshold);
    if(cuts.pt_min>0) breaks.push_back(fm/fg*selection.threshold/cuts.pt_min);
  }
  if(selection.cut==SingleLeptonCut::absolute_eta_below && std::tanh(selection.threshold)<1)
    energies.push_back(cuts.pt_min*std::cosh(selection.threshold));
  if(selection.cut!=SingleLeptonCut::none)
    for(Real energy:energies) for(Real f:{fg,fm}) breaks.push_back(2*energy/(M*f));
  const auto hard = vector_photon_pair_mass_hard_coefficient(M,mu);
  result.projected = integrate([&](Real y) {
    return hard.resolved_density(y)*B(y);
  },lower,upper,options.integration,detail::resolved_cut_breakpoints(breaks,lower,upper)).value;
  if (cuts.pt_min == 0 && std::tanh(cuts.eta_max) == 1 && selection.cut==SingleLeptonCut::none) {
    result.total = result.projected;
    return result;
  }
  result.geometry = integrate([&](Real y) {
    const Real projected = B(y), span = 1-y;
    std::vector<Real> inner;
    for(Real energy:energies) {
      inner.push_back(std::sqrt(std::max(Real{0},(2*energy/(M*fg)-y)/span)));
      inner.push_back(std::sqrt(std::max(Real{0},(1-2*energy/(M*fm))/span)));
    }
    return integrate([&](Real u) {
      const Real t = span*u*u;
      const Real density = (std::pow(1-y,2)+std::pow(1-t,2))/(y*t);
      return 2*span*u*density*(acceptance(M,y,t,fg,fm,cuts,selection)-projected);
    },std::sqrt(options.collinear_cutoff),1,options.integration,
        detail::resolved_cut_breakpoints(inner,std::sqrt(options.collinear_cutoff),1)).value;
  },lower,upper,options.integration,detail::resolved_cut_breakpoints(breaks,lower,upper)).value;
  result.total = result.projected+result.geometry;
  return result;
}
} // namespace leptonff
