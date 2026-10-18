#include "leptonff/ff_decay.hpp"
#include "leptonff/xspace_evolution.hpp"
#include "reference_config.hpp"

#include <iomanip>
#include <iostream>
#include <limits>

namespace {
using namespace leptonff;
void matrix(const MomentMatrix& a) {
  std::cout<<'[';
  for(unsigned i=0;i<3;++i) { if(i)std::cout<<',';std::cout<<'[';
    for(unsigned j=0;j<3;++j){if(j)std::cout<<',';std::cout<<a[i][j];}std::cout<<']'; }
  std::cout<<']';
}
void lattice(const LatticeMatrix& a) {
  std::cout<<'[';
  for(unsigned i=0;i<3;++i) {if(i)std::cout<<',';std::cout<<'[';
    for(unsigned j=0;j<3;++j){if(j)std::cout<<',';std::cout<<'[';
      for(std::size_t k=0;k<a[i][j].size();++k){if(k)std::cout<<',';std::cout<<a[i][j][k];}
      std::cout<<']';}std::cout<<']';}
  std::cout<<']';
}
void point(const char* label,const ReferenceInput& input,const MomentEvolutionParameters& p,
           QEDXSpaceGrid grid,bool dump) {
  const QEDXSpaceEvolution e(p,grid);
  const auto match=decay_xspace_ff_match(input.mass,input.cuts,e,input.integration);
  const auto& d=e.diagnostics();
  std::cout<<"{\"label\":\""<<label<<"\",\"bins\":"<<grid.bins
    <<",\"log_span\":"<<grid.log_span<<",\"weighted_wrap_bound\":"<<d.weighted_wrap_bound
    <<",\"weighted_tail_bound\":"<<d.weighted_tail_bound<<",\"imaginary_mass\":"<<d.imaginary_mass
    <<",\"negative_weight_mass\":"<<d.negative_weight_mass
    <<",\"maximum_steps\":"<<d.maximum_steps
    <<",\"frequency_step_shift\":"<<d.maximum_frequency_step_shift
    <<",\"weighted_step_shift_bound\":"<<d.weighted_step_shift_bound
    <<",\"born\":"<<match.born<<",\"resummed_ff_pairs\":"<<match.resummed
    <<",\"charge_flipped_pairs\":"<<match.charge_flipped_pairs
    <<",\"lattice_coefficient1\":"<<match.lattice_coefficient1
    <<",\"continuum_coefficient1\":"<<match.continuum_coefficient1
    <<",\"matched_ff_pairs\":"<<match.matched<<",\"moments\":[";
  bool first=true;
  for(Real N:{2.L,3.L,5.L,8.L}) {
    if(!first)std::cout<<',';first=false;
    const auto observed=e.moment(N),reference=evolve_qed_moment(N,p).fragmentation;
    Real error=0;
    for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
      error=std::max(error,std::abs(observed[i][j]-reference[i][j]));
    std::cout<<"{\"N\":"<<N<<",\"maximum_residual\":"<<error<<",\"matrix\":";
    matrix(observed);std::cout<<'}';
  }
  std::cout<<']';
  if(dump) for(const auto& entry:{std::make_pair("kernel0",&e.kernel0()),
      std::make_pair("kernel1",&e.kernel1()),std::make_pair("boundary1",&e.boundary1()),
      std::make_pair("propagator",&e.propagator()),std::make_pair("weights",&e.weights()),
      std::make_pair("coefficient1",&e.coefficient1()),std::make_pair("coefficient2",&e.coefficient2())}) {
    std::cout<<",\""<<entry.first<<"\":";lattice(*entry.second);
  }
  std::cout<<'}'<<std::flush;
}
}
int main(int argc,char** argv) {
  try {
    if(argc<2 || argc>3) throw std::invalid_argument("Usage: leptonff_xspace config [bins]");
    ReferenceInput input;read_reference_config(argv[1],input);
    input.qed.validate();input.cuts.validate();
    QEDXSpaceGrid grid;
    if(argc==3){std::size_t used=0;grid.bins=std::stoull(argv[2],&used);
      if(used!=std::string(argv[2]).size())throw std::invalid_argument("Invalid bins");}
    grid.validate();
    if(grid.bins>32768)throw std::invalid_argument("Refinement driver requires bins<=32768");
    MomentEvolutionParameters p;
    p.alpha_msbar_at_mu0=input.qed.alpha;p.muon_mass=input.qed.muon_mass;
    p.mu_0=input.qed.mu_0;p.mu_frag=input.qed.mu_frag;
    std::cout<<std::setprecision(18)<<"{\n\"version\":\""<<LEPTONFF_VERSION
      <<"\",\n\"compiler\":\""<<__VERSION__<<"\",\n\"real_mantissa_bits\":"<<std::numeric_limits<Real>::digits
      <<",\n\"approximation\":\"energy-weighted coupled NLO x-space; two-loop MSbar coupling; Born two-direction FF pair density; no hard H1\",\n"
      <<"\"coupling_input\":\"config.alpha is MSbar alpha(mu_0); no on-shell conversion\",\n"
      <<"\"alpha_msbar_at_mu0\":"<<p.alpha_msbar_at_mu0<<",\n\"muon_mass\":"<<p.muon_mass
      <<",\n\"mu_0\":"<<p.mu_0<<",\n\"mu_frag\":"<<p.mu_frag<<",\n\"parent_mass\":"<<input.mass
      <<",\n\"pt_min\":"<<input.cuts.pt_min<<",\n\"eta_max\":"<<input.cuts.eta_max
      <<",\n\"mass_min\":"<<input.cuts.mass_min<<",\n\"mass_max\":"<<input.cuts.mass_max
      <<",\n\"points\":[\n";
    auto coarse=grid;coarse.bins=64;point("independent_coarse",input,p,coarse,true);
    std::cout<<",\n";point("baseline",input,p,grid,false);
    grid.bins*=2;std::cout<<",\n";point("grid_x2",input,p,grid,false);
    grid.bins*=2;std::cout<<",\n";point("grid_x4",input,p,grid,false);
    grid.bins*=2;grid.log_span*=2;std::cout<<",\n";point("span_x2",input,p,grid,false);
    std::cout<<"\n]}\n";
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
