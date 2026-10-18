#include "leptonff/photonic_nll.hpp"
#include "leptonff/nlo_lepton.hpp"
#include <iomanip>
#include <cmath>
#include <iostream>

int main() {
  using namespace leptonff;
  QEDParameters q;q.mu_frag=80.385L;
  std::cout<<std::setprecision(19)<<"{\"alpha\":"<<q.alpha<<",\"mass\":"<<q.muon_mass
    <<",\"mu0\":"<<q.mu_0<<",\"mu\":"<<q.mu_frag
    <<",\"source_number\":"<<photonic_source_number_moment()<<",\"rows\":[";
  bool first=true;
  for(std::size_t bins:{1024u,4096u,65536u}) {
    LLGrid g;g.bins=bins;g.log_span=40;g.soft_log_cutoff=g.spacing()/2;
    const PhotonicNLL e(q,g);
    for(Real N:{1.L,2.L,3.L,5.L}) {
      if(!first)std::cout<<',';
      first=false;
      Real evolved=0,order2=0;
      for(std::size_t k=0;k<g.bins;++k) {
        const Real test=std::pow(e.baseline().evolution().z(k),N-1);
        evolved+=e.evolved_kernel_weights()[k]*test;order2+=e.coefficient2_weights()[k]*test;
      }
      std::cout<<"{\"bins\":"<<bins<<",\"N\":"<<N<<",\"moment\":"<<e.moment(N)
        <<",\"exact\":"<<e.exact_moment(N)<<",\"kernel\":"<<unitary_photonic_nlo_kernel().moment(N)
        <<",\"transported_kernel\":"<<evolved<<",\"coefficient2\":"<<order2
        <<",\"wrap_bound\":"<<e.diagnostics().wrap_bound<<",\"tail_bound\":"<<e.diagnostics().omitted_tail_bound
        <<",\"raw_inclusive_factor\":"<<e.diagnostics().raw_inclusive_factor<<'}';
    }
  }
  std::cout<<"]}\n";
}
