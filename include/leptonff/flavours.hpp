#pragma once
#include "leptonff/nlo_qed.hpp"
#include <string>

namespace leptonff {
struct ChargedFermion {
  std::string name;
  Real charge=-1,mass=1;
  unsigned colours=1;
  Real matching_scale=0; // zero selects the pole mass; otherwise explicit GeV
  bool pointlike=false; // explicit massless coloured boundary at matching_scale
  Real threshold() const { return matching_scale==0?mass:matching_scale; }
};
// Flat parent/daughter indices: fermion f=2*f, antifermion=2*f+1,
// photon=2*fermions.size(). Quark entries are PURE QED, not QCD evolution
// or a model for low-energy hadronic vacuum polarisation.
struct QEDFlavourScheme {
  std::vector<ChargedFermion> fermions;
  std::size_t muon_flavour=0;
  void validate() const;
  bool include_gluon=false;
  std::size_t size() const { return 2*fermions.size()+1+(include_gluon?1:0); }
  std::size_t photon() const { return 2*fermions.size(); }
  std::size_t gluon() const { return photon()+1; }
  std::array<std::size_t,3> selected() const { return {2*muon_flavour,2*muon_flavour+1,photon()}; }
  std::vector<bool> active(Real scale) const;
  static QEDFlavourScheme leptons(Real heavy_threshold_factor=1);
  static QEDFlavourScheme hadronic_pointlike(Real matching_scale=2);
  // Strict lines: name charge colours pole_mass matching_scale. '#' comments.
  static QEDFlavourScheme from_file(const std::string& path);
};
struct QEDChargeSums {
  Real squared=0,fourth=0;
  Real beta0() const { return 2*squared/3; }
  Real beta1() const { return fourth; } // da/dln(mu^2)=beta0*a^2+beta1*a^3
};
QEDChargeSums qed_charge_sums(const QEDFlavourScheme&,const std::vector<bool>& active);
Distribution flavour_leading_kernel(const QEDFlavourScheme&,const std::vector<bool>&,std::size_t parent,std::size_t observed);
Distribution flavour_next_to_leading_kernel(const QEDFlavourScheme&,const std::vector<bool>&,std::size_t parent,std::size_t observed);
// Analytically cancel the pure singlet BEFORE the finite N=1 valence moment.
Distribution flavour_valence_kernel(const QEDFlavourScheme&,const std::vector<bool>&,std::size_t fermion);
Distribution flavour_initial_condition(const QEDFlavourScheme&,const std::vector<bool>&,std::size_t parent,std::size_t observed,Real mu0);
// Non-identity part of the NLO upward TIME-LIKE matching map. The new heavy
// rows are a*K_hgamma*D_gamma below, not an identity on inactive heavy rows.
Distribution flavour_threshold_coefficient(const QEDFlavourScheme&,std::size_t heavy,
    std::size_t parent,std::size_t observed,Real matching_scale);
// One-loop OS -> MSbar conversion for leptons only, with pole masses.
// Rational completion has the same NLO expansion; no hadronic input implied.
Real leptonic_alpha_msbar(Real alpha_os,Real mu,const QEDFlavourScheme&);
Distribution qcd_leading_kernel(const QEDFlavourScheme&,const std::vector<bool>&,
    std::size_t parent,std::size_t observed);
} // namespace leptonff
