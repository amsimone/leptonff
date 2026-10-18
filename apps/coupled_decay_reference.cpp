#include "leptonff/coupled_decay.hpp"
#include "reference_config.hpp"
#include <iomanip>
#include <iostream>
#include <limits>

namespace {
using namespace leptonff;
void array(const std::vector<Real>& a) {
  std::cout<<'['; bool first=true;
  for(Real x:a){if(!first)std::cout<<',';first=false;std::cout<<x;}std::cout<<']';
}
void matrix(const std::vector<std::vector<Real>>& a) {
  std::cout<<'[';bool first=true;
  for(const auto& row:a){if(!first)std::cout<<',';first=false;array(row);}std::cout<<']';
}
void point(const char* label,const CoupledFiducialDecay& hard,
    const QEDXSpaceEvolution& e,const ReferenceInput& input,std::size_t panels,bool dump) {
  const auto r=hard.match(e,input.integration);
  const auto& d=e.diagnostics();
  std::cout<<"{\"label\":"<<std::quoted(label)<<",\"bins\":"<<e.grid().bins
    <<",\"log_span\":"<<e.grid().log_span<<",\"ff_panels\":"<<panels
    <<",\"mu_frag\":"<<e.parameters().mu_frag<<",\"coupling\":"<<std::quoted(coupling_mode_name(e.parameters().coupling))
    <<",\"a_initial\":"<<r.a_initial<<",\"a_hard\":"<<r.a_hard
    <<",\"born\":"<<r.born<<",\"ff\":"<<r.ff<<",\"charge_flipped_pairs\":"<<r.charge_flipped_pairs
    <<",\"muon_projected\":"<<r.muon_projected<<",\"muon_geometry\":"<<r.muon_geometry
    <<",\"photon_projected\":"<<r.photon_projected<<",\"photon_geometry\":"<<r.photon_geometry
    <<",\"lattice_coefficient1\":"<<r.lattice_coefficient1
    <<",\"continuum_coefficient1\":"<<r.continuum_coefficient1<<",\"coefficient2\":"<<r.coefficient2
    <<",\"overlap_correction\":"<<r.overlap_correction<<",\"raw\":"<<r.raw
    <<",\"matched\":"<<r.matched<<",\"fixed_order\":"<<r.fixed_order
    <<",\"compression_moment_residual\":"<<r.compression_moment_residual
    <<",\"compression_absolute_norm\":"<<r.compression_absolute_norm
    <<",\"weighted_wrap_bound\":"<<d.weighted_wrap_bound<<",\"weighted_tail_bound\":"<<d.weighted_tail_bound
    <<",\"weighted_step_shift_bound\":"<<d.weighted_step_shift_bound;
  if(dump) {
    for(const auto& entry:{std::make_pair("fractions",&hard.fractions()),
        std::make_pair("hard_muon",&hard.muon_projection()),std::make_pair("hard_photon",&hard.photon_projection())}) {
      std::cout<<','<<std::quoted(entry.first)<<':';array(*entry.second);
    }
    std::cout<<",\"geometry_muon\":";matrix(hard.muon_geometry());
    std::cout<<",\"geometry_photon\":";matrix(hard.photon_geometry());
    std::cout<<",\"weighted_arrays\":{";
    bool first=true;
    for(const auto& entry:{std::make_pair("same",&e.weights()[0][0]),std::make_pair("opposite",&e.weights()[0][1]),
        std::make_pair("gamma",&e.weights()[2][0]),std::make_pair("c1_same",&e.coefficient1()[0][0]),
        std::make_pair("c1_opposite",&e.coefficient1()[0][1]),std::make_pair("c1_gamma",&e.coefficient1()[2][0]),
        std::make_pair("c2_same",&e.coefficient2()[0][0])}) {
      if(!first)std::cout<<',';first=false;
      std::cout<<std::quoted(entry.first)<<':';
      const auto n=hard.muon_projection().size();
      array(std::vector<Real>(entry.second->begin(),entry.second->begin()+n));
    }
    std::cout<<'}';
  }
  std::cout<<'}'<<std::flush;
}
}
int main(int argc,char** argv) {
  try {
    if(argc<2 || argc>5) throw std::invalid_argument("Usage: leptonff_coupled config [bins [ff_panels [--scan|--dump]]]");
    ReferenceInput input;read_reference_config(argv[1],input);
    QEDXSpaceGrid grid;
    CoupledDecayOptions options;
    auto integer=[](const char* s){std::size_t used=0;auto n=std::stoull(s,&used);
      if(used!=std::string(s).size())throw std::invalid_argument("Invalid integer option");return n;};
    if(argc>=3)grid.bins=integer(argv[2]);
    if(argc>=4)options.ff_panels=integer(argv[3]);
    const std::string mode=argc==5?argv[4]:"";
    if(mode!="" && mode!="--scan" && mode!="--dump")throw std::invalid_argument("Unknown coupled driver mode");
    if(mode=="--scan" && (grid.bins>32768 || options.ff_panels>32))
      throw std::invalid_argument("Refinement scan needs bins<=32768, panels<=32");
    grid.validate(); options.validate();
    MomentEvolutionParameters p;
    p.alpha_msbar_at_mu0=input.qed.alpha;p.muon_mass=input.qed.muon_mass;
    p.mu_0=input.qed.mu_0;p.mu_frag=input.qed.mu_frag;
    std::cout<<std::setprecision(18)<<"{\"version\":\""<<LEPTONFF_VERSION
      <<"\",\"compiler\":"<<std::quoted(__VERSION__)<<",\"real_digits\":"<<std::numeric_limits<Real>::digits
      <<",\"model\":\"one_family_unpolarised_vector_at_rest_inclusive_distinct_parent_pairs\""
      <<",\"accuracy\":\"Oalpha_hard_matched_coupled_NLL_with_numerical_grid_and_geometry_controls\""
      <<",\"scheme\":\"MSbar\",\"coupling_input\":\"config.alpha is MSbar alpha(mu_0), no OS conversion\""
      <<",\"alpha_msbar_at_mu0\":"<<p.alpha_msbar_at_mu0<<",\"muon_mass\":"<<p.muon_mass
      <<",\"mu_0\":"<<p.mu_0<<",\"mu_frag\":"<<p.mu_frag<<",\"parent_mass\":"<<input.mass
      <<",\"pt_min\":"<<input.cuts.pt_min<<",\"eta_max\":"<<input.cuts.eta_max
      <<",\"mass_min\":"<<input.cuts.mass_min<<",\"mass_max\":"<<input.cuts.mass_max
      <<",\"soft_cutoff\":"<<options.geometry.soft_cutoff
      <<",\"collinear_cutoff\":"<<options.geometry.collinear_cutoff<<",\"points\":[";
    const QEDXSpaceEvolution base(p,grid);
    const CoupledFiducialDecay hard(input.mass,input.cuts,grid,options);
    point("baseline",hard,base,input,options.ff_panels,mode=="--dump");
    if(mode=="--scan") {
      options.ff_panels*=2;
      const CoupledFiducialDecay ff(input.mass,input.cuts,grid,options);
      std::cout<<',';point("ff_panels_x2",ff,base,input,options.ff_panels,false);
      grid.bins*=2;
      const QEDXSpaceEvolution fine(p,grid);
      const CoupledFiducialDecay gh(input.mass,input.cuts,grid,options);
      std::cout<<',';point("grid_x2",gh,fine,input,options.ff_panels,false);
      grid.bins*=2;grid.log_span*=2;
      const QEDXSpaceEvolution span(p,grid);
      const CoupledFiducialDecay sh(input.mass,input.cuts,grid,options);
      std::cout<<',';point("span_x2",sh,span,input,options.ff_panels,false);
      options.geometry.soft_cutoff/=4;options.geometry.collinear_cutoff/=4;
      options.geometry.integration.absolute_tolerance/=4;
      options.geometry.integration.relative_tolerance/=4;
      const CoupledFiducialDecay final(input.mass,input.cuts,grid,options);
      std::cout<<',';point("hard_slivers",final,span,input,options.ff_panels,false);
      for(Real scale:{.5L,2.0L}) {
        p.mu_frag=scale*input.qed.mu_frag;
        std::cout<<',';point(scale<1?"scale_half":"scale_double",final,QEDXSpaceEvolution(p,grid),input,options.ff_panels,false);
      }
    }
    std::cout<<"]}\n";
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
