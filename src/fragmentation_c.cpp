#include "leptonff/fragmentation_c.h"
#include "leptonff/core.hpp"
#include "api_error.hpp"
#include <cmath>
#include <memory>
#include <stdexcept>

using namespace leptonff;
using leptonff::detail::guarded;
namespace {
Species species(int id) {
  if (id<0 || id>2) throw std::invalid_argument("Invalid fragmentation species");
  return static_cast<Species>(id);
}
MomentEvolutionParameters parameters(const leptonff_fragmentation_config& c) {
  if (c.struct_size!=sizeof(c) || c.abi_version!=LEPTONFF_FF_ABI_VERSION)
    throw std::invalid_argument("Unsupported fragmentation configuration ABI");
  if (c.model!=LEPTONFF_FF_ORDER1 && c.model!=LEPTONFF_FF_COUPLED_NLL)
    throw std::invalid_argument("Unsupported fragmentation model");
  if (c.running<0 || c.running>2 || c.flavours<0 || c.flavours>1)
    throw std::invalid_argument("Unsupported running or flavour setup");
  EvolutionTheory theory;
  theory.alpha={c.alpha,static_cast<AlphaScheme>(c.alpha_scheme),c.alpha_reference_scale};
  theory.muon_mass=c.muon_mass; theory.initial_scale=c.initial_scale; theory.final_scale=c.final_scale;
  theory.running=static_cast<QEDCouplingMode>(c.running);
  if (c.flavours==LEPTONFF_LEPTON_THRESHOLDS) {
    theory.flavours=QEDFlavourScheme::leptons();
    theory.flavours->fermions[theory.flavours->muon_flavour].mass=c.muon_mass;
  }
  if (c.model==LEPTONFF_FF_ORDER1 && c.flavours!=LEPTONFF_ONE_FAMILY)
    throw std::invalid_argument("Fixed-order interface uses one-family boundaries; select coupled NLL for thresholds");
  return theory.parameters();
}
}
struct leptonff_fragmentation_context {
  int model;
  MomentEvolutionParameters theory;
  IntegrationOptions integration;
  MomentEvolutionOptions solver;
  std::unique_ptr<QEDXSpaceEvolution> lattice;
  explicit leptonff_fragmentation_context(const leptonff_fragmentation_config& c)
      :model(c.model),theory(parameters(c)) {
    if (!std::isfinite(c.absolute_tolerance) || c.absolute_tolerance<=0 ||
        !std::isfinite(c.relative_tolerance) || c.relative_tolerance<=0)
      throw std::invalid_argument("Integration tolerances must be positive and finite");
    integration.absolute_tolerance=c.absolute_tolerance;
    integration.relative_tolerance=c.relative_tolerance;
    solver.kernel_integration=integration;
    if (model==LEPTONFF_FF_COUPLED_NLL && c.grid_bins) {
      QEDXSpaceGrid grid; grid.bins=c.grid_bins; grid.log_span=c.log_span;
      lattice=std::make_unique<QEDXSpaceEvolution>(theory,grid,solver);
    }
  }
  Distribution fixed_order(Species parent,Species observed) const {
    return fragmentation_order1(parent,observed,
      {theory.alpha_msbar_at_mu0,theory.muon_mass,theory.mu_0,theory.mu_frag});
  }
};
extern "C" {
void leptonff_fragmentation_config_default(leptonff_fragmentation_config* c) {
  if (!c) return;
  *c={}; c->struct_size=sizeof(*c); c->abi_version=LEPTONFF_FF_ABI_VERSION;
  c->model=LEPTONFF_FF_ORDER1; c->alpha_scheme=LEPTONFF_ALPHA_MSBAR;
  c->running=LEPTONFF_ALPHA_TWO_LOOP; c->alpha=1/137.035999084;
  c->muon_mass=.1056583755; c->initial_scale=c->alpha_reference_scale=c->muon_mass;
  c->final_scale=91.1876; c->log_span=40;
  c->absolute_tolerance=c->relative_tolerance=1e-10;
}
int leptonff_fragmentation_create(const leptonff_fragmentation_config* c,leptonff_fragmentation_context** out) {
  if (out) *out=nullptr;
  return guarded([&] {
    if (!c || !out) throw std::invalid_argument("Null fragmentation configuration/output");
    *out=new leptonff_fragmentation_context(*c);
  });
}
void leptonff_fragmentation_destroy(leptonff_fragmentation_context* c) { delete c; }
int leptonff_fragmentation_moment(const leptonff_fragmentation_context* c,int i,int j,double N,double* out) {
  return guarded([&] {
    if (!c || !out) throw std::invalid_argument("Null fragmentation context/output");
    const auto parent=species(i),observed=species(j);
    const Real result=c->model==LEPTONFF_FF_ORDER1 ? c->fixed_order(parent,observed).moment(N,c->integration)
      : evolve_qed_moment(N,c->theory,c->solver).fragmentation[i][j];
    *out=static_cast<double>(result);
  });
}
int leptonff_fragmentation_integral(const leptonff_fragmentation_context* c,int i,int j,
                                   double lower,double upper,double* out) {
  return guarded([&] {
    if (!c || !out) throw std::invalid_argument("Null fragmentation context/output");
    const auto parent=species(i),observed=species(j);
    if (!std::isfinite(lower) || !std::isfinite(upper) || lower<0 || lower>=upper || upper>1)
      throw std::invalid_argument("Fragmentation bin requires 0<=lower<upper<=1");
    if (c->model==LEPTONFF_FF_ORDER1) *out=static_cast<double>(c->fixed_order(parent,observed).integral(lower,upper,c->integration));
    else {
      if (!c->lattice) throw std::invalid_argument("Coupled bin integrals require grid_bins>0 at creation");
      *out=static_cast<double>(c->lattice->apply(parent,observed,
        [upper](Real z) { return (z<upper || (upper==1 && z==1)) ? 1.0L : 0.0L; },lower));
    }
  });
}
int leptonff_fragmentation_convolve(const leptonff_fragmentation_context* c,int i,int j,
                                   double x,leptonff_function f,void* data,double* out) {
  return guarded([&] {
    if (!c || !out || !f) throw std::invalid_argument("Null convolution input/output");
    const auto parent=species(i),observed=species(j);
    if (c->model!=LEPTONFF_FF_ORDER1) throw std::invalid_argument("This convolution interface requires fixed order");
    *out=static_cast<double>(c->fixed_order(parent,observed).convolve(
      [=](Real z) { return static_cast<Real>(f(static_cast<double>(z),data)); },x,c->integration));
  });
}
}
