#include "leptonff/xspace_evolution.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
template<class Matrix> void matrix(const Matrix& m) {
  std::cout<<'[';bool first=true;for(const auto& row:m) {
    if(!first)std::cout<<',';first=false;std::cout<<'[';
    bool entry=true;for(auto x:row){if(!entry)std::cout<<',';entry=false;std::cout<<x;}std::cout<<']';
  }std::cout<<']';
}
}
int main(int argc,char** argv) {
  try {
    using namespace leptonff;
    if(argc>2)throw std::invalid_argument("Usage: leptonff_flavours [xspace_bins; zero skips FFT]");
    const std::size_t bins=argc==2?std::stoul(argv[1]):256;
    std::cout<<std::setprecision(18)<<"{\"version\":\""<<LEPTONFF_VERSION<<"\",\"cases\":[";
    bool first=true;
    for(bool coloured:{false,true})for(Real factor:{.5L,1.L,2.L})for(auto mode:{QEDCouplingMode::one_loop,QEDCouplingMode::two_loop}) {
      MomentEvolutionParameters p;p.flavour_scheme=QEDFlavourScheme::leptons(factor);
      if(coloured) {
        QEDFlavourScheme s;s.fermions={{"mu",-1,p.muon_mass,1,0},{"test_up",2.L/3,1.5L,3,1.5L*factor},
          {"test_down",-1.L/3,4.75L,3,4.75L/factor}};p.flavour_scheme=s;
      }else p.alpha_msbar_at_mu0=leptonic_alpha_msbar(p.alpha_msbar_at_mu0,p.mu_0,*p.flavour_scheme);
      p.coupling=mode;const auto& s=*p.flavour_scheme;const Real N=3;
      const auto r=evolve_qed_moment(N,p);
      if(!first)std::cout<<',';first=false;
      std::cout<<"{\"profile\":\""<<(coloured?"charged_test":"leptons")<<"\",\"N\":"<<N
        <<",\"coupling\":\""<<coupling_mode_name(mode)<<"\",\"alpha_msbar_at_mu0\":"<<p.alpha_msbar_at_mu0
        <<",\"mu0\":"<<p.mu_0<<",\"mu_frag\":"<<p.mu_frag<<",\"muon_index\":"<<s.muon_flavour<<",\"fermions\":[";
      bool ffirst=true;for(const auto& f:s.fermions) {
        if(!ffirst)std::cout<<',';ffirst=false;
        std::cout<<"{\"name\":\""<<f.name<<"\",\"charge\":"<<f.charge<<",\"colours\":"<<f.colours
          <<",\"mass\":"<<f.mass<<",\"threshold\":"<<f.threshold()<<'}';
      }
      std::cout<<"],\"a_final\":"<<r.a_final<<",\"beta0_integral\":"<<run_flavour_coupling(p,p.mu_frag).beta0_integral
        <<",\"steps_per_interval\":"<<r.steps<<",\"step_shift\":"<<r.step_refinement_shift;
      for(const auto& entry:{std::make_pair("fragmentation",r.flavour_fragmentation),std::make_pair("evolution",r.flavour_evolution),
                            std::make_pair("coefficient1",r.flavour_coefficient1),std::make_pair("coefficient2",r.flavour_coefficient2)}) {
        std::cout<<",\""<<entry.first<<"\":";matrix(entry.second);
      }
      const auto active=s.active(p.mu_frag);std::vector<std::vector<Real>> k0(s.size(),std::vector<Real>(s.size())),k1=k0;
      for(std::size_t i=0;i<s.size();++i)for(std::size_t j=0;j<s.size();++j) {
        k0[i][j]=flavour_leading_kernel(s,active,i,j).moment(N);
        k1[i][j]=flavour_next_to_leading_kernel(s,active,i,j).moment(N);
      }
      std::cout<<",\"final_kernel0\":";matrix(k0);std::cout<<",\"final_kernel1\":";matrix(k1);
      if(bins && factor==1 && mode==QEDCouplingMode::two_loop) {
        QEDXSpaceGrid g;g.bins=bins;const QEDXSpaceEvolution x(p,g);
        std::cout<<",\"lattice_moment2\":";matrix(x.flavour_moment(2));
        std::cout<<",\"continuum_moment2\":";matrix(evolve_qed_moment(2,p).flavour_fragmentation);
        std::cout<<",\"lattice_moment3\":";matrix(x.flavour_moment(3));
        std::cout<<",\"lattice_bins\":"<<bins<<",\"wrap_bound\":"<<x.diagnostics().weighted_wrap_bound
          <<",\"tail_bound\":"<<x.diagnostics().weighted_tail_bound;
      }
      std::cout<<'}';
    }
    std::cout<<"]}\n";
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
