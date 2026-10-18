#include "leptonff/c_api.h"
#include "api_error.hpp"
#include "leptonff/production_decay.hpp"
#include "leptonff/hadronic_input.hpp"
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <string>

using namespace leptonff;
namespace {
using leptonff::detail::guarded;
bool neutral(int p) { return p==LEPTONFF_PROCESS_Z || p==LEPTONFF_Z; }
bool positive(int p) { return p==LEPTONFF_PROCESS_W_PLUS || p==LEPTONFF_WPLUS; }
FourVector vector(const double* p) {if(!p)throw std::invalid_argument("Null momentum");return {p[0],p[1],p[2],p[3]};}
DecayConfiguration configuration(const leptonff_decay_config& c) {
  if(c.struct_size!=sizeof(c) || c.abi_version!=LEPTONFF_ABI_VERSION)throw std::invalid_argument("Unsupported leptonff configuration ABI");
  DecayConfiguration q;
  if(!positive(c.process) && !neutral(c.process) && c.process!=LEPTONFF_PROCESS_W_MINUS && c.process!=LEPTONFF_WMINUS)throw std::invalid_argument("Unsupported W/Z process");
  q.process=neutral(c.process)?SpinDecay::neutral_vector:positive(c.process)?SpinDecay::wplus:SpinDecay::wminus;
  q.evolution.muon_mass=c.muon_mass;q.evolution.mu_0=c.muon_mass;
  q.evolution.mu_frag=c.mu_frag;q.evolution.alpha_msbar_at_mu0=c.alpha;
  q.grid.bins=c.grid_bins;q.grid.log_span=c.log_span;q.hard_samples=c.hard_samples;
  q.soft_cutoff=c.soft_cutoff;q.collinear_cutoff=c.collinear_cutoff;q.minimum_z=c.minimum_z;
  if(c.model==LEPTONFF_PHOTONIC_LL || c.model==LEPTONFF_PHOTONIC_NLL)
    q.mode=c.model==LEPTONFF_PHOTONIC_LL?DecayMode::photonic_ll:DecayMode::photonic_nll;
  else if(c.model==LEPTONFF_LEPTONIC_NLL) {
    q.mode=DecayMode::coupled_nll;q.evolution.flavour_scheme=QEDFlavourScheme::leptons();
    q.evolution.alpha_msbar_at_mu0=leptonic_alpha_msbar(c.alpha,c.muon_mass,*q.evolution.flavour_scheme);
    q.fraction_measure=FractionMeasure::conservative_cells;
  }else if(c.model==LEPTONFF_QED_HVP) {
    q.mode=DecayMode::coupled_nll;configure_hadronic_pointlike(q.evolution,HadronicAlphaInput{},c.alpha,c.hadronic_matching_scale);
    q.evolution.qcd->alpha_s_reference=0;q.fraction_measure=FractionMeasure::conservative_cells;
  }else throw std::invalid_argument("Unsupported decay model");
  return q;
}
LabSpectrum spectrum(const leptonff_decay_config& c,const leptonff_observable* obs,size_t count) {
  WLabCuts w;w.muon_pt_min=c.muon_pt_min;w.muon_eta_max=c.muon_eta_max;w.neutrino_pt_min=c.neutrino_pt_min;w.transverse_mass_min=c.transverse_mass_min;
  FiducialCuts z;z.pt_min=c.muon_pt_min;z.eta_max=c.muon_eta_max;z.mass_min=c.mass_min;z.mass_max=c.mass_max;
  const bool is_neutral=neutral(c.process);
  auto bins=LabSpectrum::reference_binning(is_neutral);
  if(count) {
    if(!obs)throw std::invalid_argument("Null observable definitions");
    bins.clear();
    for(size_t k=0;k<count;++k) {
      if(!obs[k].observable || !obs[k].edges || obs[k].edge_count<2)throw std::invalid_argument("Invalid observable definition");
      bins.push_back({obs[k].observable,{}});
      for(size_t i=0;i<obs[k].edge_count;++i)bins.back().edges.push_back(obs[k].edges[i]);
    }
  }
  return LabSpectrum::configured(is_neutral,bins,w,z);
}

}
struct leptonff_decay_context {
  DecayCalculator calculator;LabSpectrum selections;
  leptonff_decay_context(const leptonff_decay_config& c,const leptonff_observable* o,size_t n)
    :calculator(configuration(c)),selections(spectrum(c,o,n)) {}
};
extern "C" {
void leptonff_decay_config_default(leptonff_decay_config* c) {
  if(!c)return;
  *c={};c->struct_size=sizeof(*c);c->abi_version=LEPTONFF_ABI_VERSION;c->process=LEPTONFF_WPLUS;c->model=LEPTONFF_PHOTONIC_NLL;
  c->alpha=1/137.035999084;c->muon_mass=.1056583755;c->mu_frag=80.385;c->grid_bins=4096;c->log_span=40;c->hard_samples=16;
  c->soft_cutoff=1e-5;c->collinear_cutoff=1e-5;c->minimum_z=1e-4;c->muon_pt_min=25;c->muon_eta_max=2.4;
  c->mass_min=60;c->mass_max=120;c->hadronic_matching_scale=2;
}
int leptonff_decay_create(const leptonff_decay_config* c,const leptonff_observable* o,size_t n,leptonff_decay_context** out) {
  if(out)*out=nullptr;
  return guarded([&]{if(!c || !out)throw std::invalid_argument("Null context configuration/output");*out=new leptonff_decay_context(*c,o,n);});
}
void leptonff_decay_destroy(leptonff_decay_context* c) {delete c;}
size_t leptonff_decay_output_count(const leptonff_decay_context* c) {return c?1+c->selections.bins.size():0;}
int leptonff_decay_get_diagnostics(const leptonff_decay_context* c,leptonff_decay_diagnostics* out) {
  return guarded([&]{
    if(!c || !out)throw std::invalid_argument("Null diagnostics input/output");
    *out={};out->abi_version=LEPTONFF_ABI_VERSION;out->a_initial=double(c->calculator.fragmentation().a_initial());out->raw_inclusive_factor=1;
    const auto* p=c->calculator.photonic();
    if(p && p->nll_evolution()) {
      const auto& d=p->nll_evolution()->diagnostics();out->photonic_unitary_conditioning=1;
      out->source_number_moment=double(d.source_number_moment);out->raw_inclusive_factor=double(d.raw_inclusive_factor);
      out->wrap_bound=double(d.wrap_bound);out->tail_bound=double(d.omitted_tail_bound);
    }else if(c->calculator.coupled()) {
      const auto& d=c->calculator.coupled()->evolution().diagnostics();out->wrap_bound=double(d.weighted_wrap_bound);out->tail_bound=double(d.weighted_tail_bound);
    }
  });
}
int leptonff_decay_evaluate(const leptonff_decay_context* c,const double* a,const double* b,uint64_t seed,leptonff_prediction* out,size_t capacity) {
  return guarded([&]{
    if(!c || !out || capacity<leptonff_decay_output_count(c))throw std::invalid_argument("Insufficient prediction output capacity");
    const auto value=c->calculator.evaluate({vector(a),vector(b),seed},c->selections);
    const auto store=[&](size_t i,const DecayPrediction& p){out[i]={double(p.born),double(p.order1),double(p.matched),double(p.order2),double(p.higher)};};
    store(0,value.total);for(size_t k=0;k<value.bins.size();++k)store(k+1,value.bins[k]);
  });
}
int leptonff_decay_response(const leptonff_decay_context* c,const double* q,uint64_t seed,unsigned ncos,unsigned nphi,double* out,size_t capacity) {
  return guarded([&]{
    if(!c || !out || capacity<90*leptonff_decay_output_count(c))throw std::invalid_argument("Insufficient tensor output capacity");
    const auto response=c->calculator.response(vector(q),c->selections,seed,{ncos,nphi});size_t offset=0;
    const auto store=[&](const std::array<SpinTensor,5>& tensors){for(const auto& t:tensors)for(const auto& row:t)for(auto v:row){out[offset++]=double(v.real());out[offset++]=double(v.imag());}};
    store(response.total);for(const auto& b:response.bins)store(b);
  });
}
}
