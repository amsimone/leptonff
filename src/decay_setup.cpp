#include "leptonff/decay_setup.hpp"
#include "leptonff/hadronic_input.hpp"
#include "leptonff/version.hpp"
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace leptonff {
DecayPreset decay_preset(const std::string& name) {
  if(name=="born")return DecayPreset::born;
  if(name=="photonic_ll")return DecayPreset::photonic_ll;
  if(name=="photonic_nll")return DecayPreset::photonic_nll;
  if(name=="leptonic_nll")return DecayPreset::leptonic_nll;
  if(name=="qed_hvp")return DecayPreset::qed_hvp;
  throw std::invalid_argument("Unknown decay preset: "+name);
}
DecayConfiguration make_decay_configuration(SpinDecay process,DecayPreset model,const AlphaInput& alpha,
    Real mass,Real scale,Real hadronic_scale) {
  alpha.validate();DecayConfiguration c;c.process=process;
  EvolutionTheory theory;theory.alpha=alpha;theory.muon_mass=mass;theory.initial_scale=mass;theory.final_scale=scale;
  if(model==DecayPreset::born || model==DecayPreset::photonic_ll || model==DecayPreset::photonic_nll) {
    // A fixed input is a number; no OS/MSbar conversion is part of this sector.
    c.mode=model==DecayPreset::photonic_ll?DecayMode::photonic_ll:DecayMode::photonic_nll;
    c.evolution.muon_mass=mass;c.evolution.mu_0=mass;c.evolution.mu_frag=scale;
    c.evolution.alpha_msbar_at_mu0=model==DecayPreset::born?0:alpha.value;
    c.evolution.coupling=QEDCouplingMode::fixed;
  }else if(model==DecayPreset::leptonic_nll) {
    theory.flavours=QEDFlavourScheme::leptons();c.evolution=theory.parameters();
    c.mode=DecayMode::coupled_nll;c.fraction_measure=FractionMeasure::conservative_cells;
  }else if(model==DecayPreset::qed_hvp) {
    if(alpha.scheme!=AlphaScheme::on_shell)throw std::invalid_argument("HVP preset requires its on-shell alpha anchor");
    c.evolution.muon_mass=mass;c.evolution.mu_0=mass;c.evolution.mu_frag=scale;
    configure_hadronic_pointlike(c.evolution,HadronicAlphaInput{},alpha.value,hadronic_scale);
    c.evolution.qcd->alpha_s_reference=0;c.grid.log_span=64;
    c.mode=DecayMode::coupled_nll;c.fraction_measure=FractionMeasure::conservative_cells;
  }else throw std::invalid_argument("Invalid decay preset");
  c.validate();return c;
}
std::string decay_configuration_identity(const DecayConfiguration& c) {
  c.validate();const auto& p=c.evolution;std::ostringstream s;s.imbue(std::locale::classic());
  s<<std::setprecision(std::numeric_limits<Real>::max_digits10)<<"leptonff-decay-1 "<<version()<<' '
    <<std::numeric_limits<Real>::digits<<' '<<int(c.process)<<' '<<int(c.mode)<<' '
    <<p.alpha_msbar_at_mu0<<' '<<p.muon_mass<<' '<<p.mu_0<<' '<<p.mu_frag<<' '<<int(p.coupling)<<' '
    <<bool(p.flavour_scheme)<<' ';
  if(p.flavour_scheme) {
    const auto& f=*p.flavour_scheme;s<<f.muon_flavour<<' '<<f.include_gluon<<' '<<f.fermions.size()<<' ';
    for(const auto& v:f.fermions)s<<std::quoted(v.name)<<' '<<v.charge<<' '<<v.mass<<' '<<v.colours<<' '<<v.matching_scale<<' '<<v.pointlike<<' ';
  }
  s<<bool(p.qcd)<<' ';if(p.qcd)s<<p.qcd->alpha_s_reference<<' '<<p.qcd->reference_scale<<' '<<p.qcd->start_scale<<' ';
  s<<bool(p.hadronic_matching)<<' ';if(p.hadronic_matching)s<<p.hadronic_matching->scale<<' '<<p.hadronic_matching->coefficient<<' ';
  s<<c.grid.bins<<' '<<c.grid.log_span<<' '<<c.grid.wrap_exponent<<' '<<c.grid.max_wrap_bound<<' '<<c.grid.max_tail_bound<<' '
    <<c.hard_samples<<' '<<c.soft_cutoff<<' '<<c.collinear_cutoff<<' '<<c.minimum_z<<' '<<c.dress_photonic_hard<<' '<<int(c.fraction_measure);
  return s.str();
}
std::string decay_configuration_json(const DecayConfiguration& c) {
  const auto& p=c.evolution;std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(std::numeric_limits<Real>::max_digits10);
  const bool photonic=c.mode==DecayMode::photonic_ll || c.mode==DecayMode::photonic_nll;
  const int process=c.process==SpinDecay::neutral_vector?23:c.process==SpinDecay::wplus?24:-24;
  const char* mode=c.mode==DecayMode::photonic_ll?"photonic_ll":c.mode==DecayMode::photonic_nll?"photonic_nll":c.mode==DecayMode::coupled_nll?"coupled_nll":"ff_only";
  // The identity is opaque. Quote arbitrary flavour names safely for JSON.
  s<<"{\"version\":\""<<version()<<"\",\"alpha_msbar_mu0\":"<<p.alpha_msbar_at_mu0
    <<",\"mu0\":"<<p.mu_0<<",\"muon_mass\":"<<p.muon_mass<<",\"mu_frag\":"<<p.mu_frag
    <<",\"coupling\":\""<<(photonic?"fixed":coupling_mode_name(p.coupling))<<"\",\"process\":"<<process<<",\"mode\":\""<<mode<<'"'
    <<",\"resolved_alpha_scheme\":\""<<(photonic?"fixed":"msbar")<<'"'
    <<",\"grid_bins\":"<<c.grid.bins<<",\"log_span\":"<<c.grid.log_span<<",\"minimum_z\":"<<c.minimum_z
    <<",\"hard_samples\":"<<c.hard_samples<<",\"identity\":\"";
  for(unsigned char ch:decay_configuration_identity(c)) {
    if(ch=='"' || ch=='\\')s<<'\\'<<char(ch);
    else if(ch<32)s<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<unsigned(ch)<<std::dec;
    else s<<char(ch);
  }
  s<<"\"}";return s.str();
}
} // namespace leptonff
