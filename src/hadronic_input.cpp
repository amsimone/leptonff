#include "leptonff/hadronic_input.hpp"
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

namespace leptonff {
void HadronicAlphaInput::validate() const {
  if(!std::isfinite(reference_scale) || reference_scale<50 || reference_scale>120 ||
     !std::isfinite(delta_had5_os) || delta_had5_os<=0 || delta_had5_os>.1L ||
     !std::isfinite(uncertainty) || uncertainty<0 || uncertainty>=delta_had5_os ||
     !std::isfinite(alpha_s_reference) || alpha_s_reference<0 || alpha_s_reference>.3L)
    throw std::invalid_argument("Invalid five-flavour high-scale on-shell HVP anchor");
}
HadronicAlphaInput HadronicAlphaInput::from_file(const std::string& path) {
  std::ifstream f(path);if(!f)throw std::invalid_argument("Cannot read hadronic input");
  HadronicAlphaInput x;std::map<std::string,Real*> fields{{"reference_scale",&x.reference_scale},
    {"delta_had5_os",&x.delta_had5_os},{"uncertainty",&x.uncertainty},{"alpha_s_reference",&x.alpha_s_reference}};
  std::string line;std::size_t count=0;
  while(std::getline(f,line)) {
    line=line.substr(0,line.find('#'));std::istringstream row(line);std::string key,extra;Real value;
    if(!(row>>key))continue;
    const auto it=fields.find(key);if(it==fields.end() || !(row>>value) || row>>extra)throw std::invalid_argument("Invalid/duplicate hadronic input field");
    *it->second=value;fields.erase(it);++count;
  }
  if(count!=4)throw std::invalid_argument("Hadronic input requires all four fields");
  x.validate();return x;
}
Real hadronic_alpha_msbar_reference(Real alpha,const HadronicAlphaInput& x) {
  x.validate();if(!std::isfinite(alpha) || alpha<0 || alpha>=1)throw std::invalid_argument("Invalid on-shell alpha");
  if(alpha==0)return 0;
  const Real pi=std::acos(-1.L),reference=1/137.035999084L;
  // The experimental Delta is quoted at the reference alpha(0). Scale the
  // leading vacuum-polarisation contribution when testing alpha -> 0.
  MomentEvolutionParameters leptons;
  leptons.flavour_scheme=QEDFlavourScheme::leptons();
  leptons.alpha_msbar_at_mu0=leptonic_alpha_msbar(alpha,leptons.mu_0,*leptons.flavour_scheme);
  const Real alpha_leptonic=2*pi*run_flavour_coupling(leptons,x.reference_scale).a_final;
  Real shift=x.delta_had5_os*alpha/reference;
  shift+=alpha/pi*(55.L/27+x.alpha_s_reference/pi*(605.L/108-44.L*1.2020569031595942854L/9));
  shift+=1-alpha/alpha_leptonic;
  if(shift>=1)throw std::invalid_argument("HVP input outside perturbative range");
  return alpha/(1-shift);
}
void configure_hadronic_pointlike(MomentEvolutionParameters& p,const HadronicAlphaInput& x,Real alpha,Real scale) {
  x.validate();
  if(p.mu_0>=scale)throw std::invalid_argument("Hadronic calibration requires an initial scale below the matching scale");
  p.flavour_scheme=QEDFlavourScheme::hadronic_pointlike(scale);
  p.qcd=QCDEvolutionParameters{x.alpha_s_reference,x.reference_scale,scale};
  p.hadronic_matching.reset();
  const Real physical=alpha==0?1/137.035999084L:alpha;
  p.alpha_msbar_at_mu0=leptonic_alpha_msbar(physical,p.mu_0,QEDFlavourScheme::leptons());
  const Real below_extra=run_flavour_coupling(p,scale).a_final;
  auto high=p;high.mu_0=x.reference_scale;high.alpha_msbar_at_mu0=hadronic_alpha_msbar_reference(physical,x);
  const Real target=run_flavour_coupling(high,scale).a_final;
  p.hadronic_matching=HadronicMatching{scale,1/below_extra-1/target};
  if(alpha==0)p.alpha_msbar_at_mu0=0;
  p.validate();
}
} // namespace leptonff
