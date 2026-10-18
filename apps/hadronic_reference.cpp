#include "leptonff/hadronic_input.hpp"
#include "leptonff/xspace_evolution.hpp"
#include <iostream>
#include <iomanip>
namespace {
template<class Matrix> void matrix(const Matrix& a) {
  std::cout<<'[';bool first=true;for(const auto& row:a){if(!first)std::cout<<',';first=false;std::cout<<'[';bool col=true;
    for(const auto& x:row){if(!col)std::cout<<',';col=false;std::cout<<x;}std::cout<<']';}std::cout<<']';
}
}
int main(int argc,char** argv) {
  try {
    using namespace leptonff;const std::size_t bins=argc>1?std::stoul(argv[1]):128;
    const HadronicAlphaInput input;std::cout<<std::setprecision(18)<<"{\"cases\":[";bool first=true;
    for(Real N:{2.L,3.L})for(bool strong:{false,true}) {
      MomentEvolutionParameters p;configure_hadronic_pointlike(p,input,1/137.035999084L);
      if(!strong)p.qcd->alpha_s_reference=0;
      const auto& s=*p.flavour_scheme;const auto r=evolve_qed_moment(N,p);
      if(!first)std::cout<<',';first=false;
      std::cout<<"{\"N\":"<<N<<",\"mu0\":"<<p.mu_0<<",\"mu_frag\":"<<p.mu_frag
        <<",\"alpha_msbar_at_mu0\":"<<p.alpha_msbar_at_mu0<<",\"a_final\":"<<r.a_final
        <<",\"alpha_s_reference\":"<<p.qcd->alpha_s_reference<<",\"reference_scale\":"<<p.qcd->reference_scale
        <<",\"start_scale\":"<<p.qcd->start_scale<<",\"hvp_scale\":"<<p.hadronic_matching->scale
        <<",\"hvp_coefficient\":"<<p.hadronic_matching->coefficient<<",\"muon_index\":"<<s.muon_flavour<<",\"fermions\":[";
      bool f=true;for(const auto& q:s.fermions){if(!f)std::cout<<',';f=false;
        std::cout<<"{\"charge\":"<<q.charge<<",\"colours\":"<<q.colours<<",\"mass\":"<<q.mass<<",\"threshold\":"<<q.threshold()<<",\"pointlike\":"<<(q.pointlike?"true":"false")<<'}';}
      std::cout<<"],\"fragmentation\":";matrix(r.flavour_fragmentation);std::cout<<",\"evolution\":";matrix(r.flavour_evolution);
      if(bins && strong){QEDXSpaceGrid g;g.bins=bins;g.log_span=64;const QEDXSpaceEvolution x(p,g);std::cout<<",\"lattice_moment\":";matrix(x.flavour_moment(N));}
      std::cout<<'}';
    }
    std::cout<<"]}\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
