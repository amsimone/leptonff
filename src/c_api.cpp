#include "leptonff/c_api.h"
#include "api_error.hpp"
#include "leptonff/production_decay.hpp"
#include "leptonff/decay_setup.hpp"
#include <cstring>
#include <sstream>
#include <limits>
#include <iomanip>
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
DecayConfiguration configuration(const leptonff_decay_config& c,const leptonff_decay_theory* t=nullptr) {
  if(c.struct_size!=sizeof(c) || c.abi_version!=LEPTONFF_ABI_VERSION)throw std::invalid_argument("Unsupported leptonff configuration ABI");
  if(!positive(c.process) && !neutral(c.process) && c.process!=LEPTONFF_PROCESS_W_MINUS && c.process!=LEPTONFF_WMINUS)throw std::invalid_argument("Unsupported W/Z process");
  const auto process=neutral(c.process)?SpinDecay::neutral_vector:positive(c.process)?SpinDecay::wplus:SpinDecay::wminus;
  DecayPreset model;
  if(c.model==LEPTONFF_PHOTONIC_LL)model=DecayPreset::photonic_ll;
  else if(c.model==LEPTONFF_PHOTONIC_NLL)model=DecayPreset::photonic_nll;
  else if(c.model==LEPTONFF_LEPTONIC_NLL)model=DecayPreset::leptonic_nll;
  else if(c.model==LEPTONFF_QED_HVP)model=DecayPreset::qed_hvp;
  else throw std::invalid_argument("Unsupported decay model");
  AlphaInput alpha=AlphaInput::on_shell(c.alpha);
  if(t) {
    if(t->struct_size!=sizeof(*t) || t->abi_version!=1 || t->alpha_scheme<0 || t->alpha_scheme>1 || t->running<0 || t->running>2)
      throw std::invalid_argument("Invalid explicit decay theory");
    alpha={t->alpha,static_cast<AlphaScheme>(t->alpha_scheme),t->alpha_reference_scale};
  }
  const auto preset_alpha=t && model==DecayPreset::leptonic_nll?AlphaInput::msbar(alpha.value,c.muon_mass):alpha;
  auto q=make_decay_configuration(process,model,preset_alpha,c.muon_mass,c.mu_frag,c.hadronic_matching_scale);
  if(t) {
    if(q.mode!=DecayMode::coupled_nll) {
      if(t->running!=0 || t->initial_scale!=c.muon_mass || (alpha.scheme==AlphaScheme::msbar && alpha.reference_scale!=t->initial_scale))
        throw std::invalid_argument("Photonic theory requires fixed alpha and initial_scale=muon_mass");
      q.evolution.coupling=QEDCouplingMode::fixed;
    }else if(model==DecayPreset::qed_hvp) {
      if(t->initial_scale!=c.muon_mass)throw std::invalid_argument("HVP anchor requires initial_scale=muon_mass");
      q.evolution.coupling=static_cast<QEDCouplingMode>(t->running);
    }else {
      EvolutionTheory theory;theory.alpha=alpha;theory.muon_mass=c.muon_mass;theory.initial_scale=t->initial_scale;
      theory.final_scale=c.mu_frag;theory.running=static_cast<QEDCouplingMode>(t->running);theory.flavours=q.evolution.flavour_scheme;
      q.evolution=theory.parameters();
    }
  }
  q.grid.bins=c.grid_bins;q.grid.log_span=c.log_span;q.hard_samples=c.hard_samples;
  q.soft_cutoff=c.soft_cutoff;q.collinear_cutoff=c.collinear_cutoff;q.minimum_z=c.minimum_z;
  q.validate();return q;
}
void string_result(const std::string& value,char* out,size_t* capacity) {
  if(!capacity)throw std::invalid_argument("Null string capacity");
  const auto available=*capacity;*capacity=value.size()+1;
  if(out) {if(available<*capacity)throw std::invalid_argument("Insufficient string capacity");std::memcpy(out,value.c_str(),*capacity);}
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
  leptonff_decay_config config;
  DecayCalculator calculator;MeasurementPlan measurement;
  leptonff_decay_context(const leptonff_decay_config& c,const leptonff_observable* o,size_t n,const leptonff_decay_theory* t=nullptr)
    :config(c),calculator(configuration(c,t)),measurement(neutral(c.process),spectrum(c,o,n)) {}
  leptonff_decay_context(const leptonff_decay_context& source,const leptonff_decay_config& c,const leptonff_observable* o,size_t n)
    :config(c),calculator(source.calculator),measurement(neutral(c.process),spectrum(c,o,n)) {}
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
void leptonff_decay_theory_default(leptonff_decay_theory* t) {
  if(!t)return;
  *t={};t->struct_size=sizeof(*t);t->abi_version=1;t->alpha_scheme=1;t->running=0;
  t->alpha=1/137.035999084;t->initial_scale=.1056583755;t->alpha_reference_scale=t->initial_scale;
}
int leptonff_decay_create_with_theory(const leptonff_decay_config* c,const leptonff_decay_theory* t,
    const leptonff_observable* o,size_t n,leptonff_decay_context** out) {
  if(out)*out=nullptr;
  return guarded([&]{if(!c || !t || !out)throw std::invalid_argument("Null theory configuration/output");*out=new leptonff_decay_context(*c,o,n,t);});
}
int leptonff_decay_with_observables(const leptonff_decay_context* source,const leptonff_decay_config* selection,
    const leptonff_observable* o,size_t n,leptonff_decay_context** out) {
  if(out)*out=nullptr;
  return guarded([&]{
    if(!source || !out)throw std::invalid_argument("Null measurement view source/output");
    const auto& c=selection?*selection:source->config;
    auto expected=source->config;
    expected.muon_pt_min=c.muon_pt_min;expected.muon_eta_max=c.muon_eta_max;
    expected.neutrino_pt_min=c.neutrino_pt_min;expected.transverse_mass_min=c.transverse_mass_min;
    expected.mass_min=c.mass_min;expected.mass_max=c.mass_max;
    // Compare named fields rather than struct padding. Explicit-theory
    // contexts keep their resolved calculator without re-resolving presets.
    if(c.struct_size!=expected.struct_size || c.abi_version!=expected.abi_version || c.process!=expected.process || c.model!=expected.model ||
       c.alpha!=expected.alpha || c.muon_mass!=expected.muon_mass || c.mu_frag!=expected.mu_frag || c.grid_bins!=expected.grid_bins ||
       c.log_span!=expected.log_span || c.hard_samples!=expected.hard_samples || c.soft_cutoff!=expected.soft_cutoff ||
       c.collinear_cutoff!=expected.collinear_cutoff || c.minimum_z!=expected.minimum_z || c.hadronic_matching_scale!=expected.hadronic_matching_scale)
      throw std::invalid_argument("A measurement view may only change cuts and binning");
    *out=new leptonff_decay_context(*source,c,o,n);
  });
}
int leptonff_decay_configuration_json(const leptonff_decay_context* c,char* out,size_t* capacity) {
  return guarded([&]{if(!c)throw std::invalid_argument("Null context");string_result(decay_configuration_json(c->calculator.configuration()),out,capacity);});
}
int leptonff_decay_output_label(const leptonff_decay_context* c,size_t index,char* out,size_t* capacity) {
  return guarded([&]{
    if(!c || index>=leptonff_decay_output_count(c))throw std::invalid_argument("Invalid output index");
    std::ostringstream s;s<<std::setprecision(std::numeric_limits<double>::max_digits10);
    if(!index)s<<"total";else {const auto& b=c->measurement.spectrum().bins[index-1];s<<b.observable<<':'<<b.lower<<':'<<b.upper;}
    string_result(s.str(),out,capacity);
  });
}
int leptonff_decay_evaluate_batch(const leptonff_decay_context* c,const double* first,const double* second,
    const uint64_t* seeds,size_t count,unsigned threads,leptonff_prediction* out,size_t capacity) {
  return guarded([&]{
    if(!c || !threads || threads>1024)throw std::invalid_argument("Invalid context/worker count");
    const auto outputs=leptonff_decay_output_count(c);
    if(count>std::numeric_limits<size_t>::max()/4 || count>std::numeric_limits<size_t>::max()/outputs || capacity<count*outputs || (count && (!first || !second || !seeds || !out)))
      throw std::invalid_argument("Invalid batch inputs/output capacity");
    std::vector<BornDecayPoint> points;points.reserve(count);
    for(size_t i=0;i<count;++i)points.push_back({vector(first+4*i),vector(second+4*i),seeds[i]});
    const auto result=c->calculator.evaluate_batch(points,c->measurement,threads);
    const auto store=[&](size_t i,const DecayPrediction& p){out[i]={double(p.born),double(p.order1),double(p.matched),double(p.order2),double(p.higher)};};
    for(size_t i=0;i<count;++i) {store(i*outputs,result[i].total);for(size_t k=0;k<result[i].bins.size();++k)store(i*outputs+k+1,result[i].bins[k]);}
  });
}
void leptonff_decay_destroy(leptonff_decay_context* c) {delete c;}
size_t leptonff_decay_output_count(const leptonff_decay_context* c) {return c?1+c->measurement.spectrum().bins.size():0;}
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
    const auto value=c->calculator.evaluate({vector(a),vector(b),seed},c->measurement);
    const auto store=[&](size_t i,const DecayPrediction& p){out[i]={double(p.born),double(p.order1),double(p.matched),double(p.order2),double(p.higher)};};
    store(0,value.total);for(size_t k=0;k<value.bins.size();++k)store(k+1,value.bins[k]);
  });
}
int leptonff_decay_response(const leptonff_decay_context* c,const double* q,uint64_t seed,unsigned ncos,unsigned nphi,double* out,size_t capacity) {
  return guarded([&]{
    if(!c || !out || capacity<90*leptonff_decay_output_count(c))throw std::invalid_argument("Insufficient tensor output capacity");
    const auto response=c->calculator.response(vector(q),c->measurement,seed,{ncos,nphi});size_t offset=0;
    const auto store=[&](const std::array<SpinTensor,5>& tensors){for(const auto& t:tensors)for(const auto& row:t)for(auto v:row){out[offset++]=double(v.real());out[offset++]=double(v.imag());}};
    store(response.total);for(const auto& b:response.bins)store(b);
  });
}
}
