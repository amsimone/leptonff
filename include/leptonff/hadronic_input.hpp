#pragma once
#include "leptonff/moment_evolution.hpp"

namespace leptonff {
// Five-flavour timelike on-shell HVP anchor, not an MSbar coupling or an
// exclusive real-pair distribution. Default: PDG 2025 review, eq. 10.20/table 10.1.
struct HadronicAlphaInput {
  Real reference_scale=91.1876L,delta_had5_os=.02760L,uncertainty=.00010L;
  Real alpha_s_reference=.1193L;
  void validate() const;
  static HadronicAlphaInput from_file(const std::string&);
};
// Fermionic five-quark + three-lepton QED coupling. Two-loop leptonic running with the existing one-loop
// low-scale input, plus massless hadronic OS/MSbar conversion through alpha*alpha_s.
// W/top loops, higher QED, higher QCD conversion and mass power terms omitted.
Real hadronic_alpha_msbar_reference(Real alpha_os,const HadronicAlphaInput&);
// Fits an EXTRA decoupling coefficient at matching_scale without counting
// the perturbative quark running twice. Requires the explicit pointlike scheme.
void configure_hadronic_pointlike(MomentEvolutionParameters&,const HadronicAlphaInput&,
    Real alpha_os,Real matching_scale=2);
} // namespace leptonff
