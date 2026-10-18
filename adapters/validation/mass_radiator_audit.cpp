// SPDX-License-Identifier: MIT
// Numerical outputs only: the independent reference algebra lives in Python.
#include "../dyturbo/mass_radiator.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc,char** argv) {try {
  using namespace leptonff;
  const unsigned cells=argc>1?std::stoul(argv[1]):65536;
  const bool profile=argc>2 && std::string(argv[2])=="profile";
  LLGrid grid;grid.bins=cells;grid.log_span=32;grid.soft_log_cutoff=grid.spacing()/2;
  std::cout<<std::setprecision(18);
  std::cout<<"{\"schema\":\"leptonff-complete-mass-audit-v1\",\"cells\":"<<cells<<",\"rows\":[";
  bool first=true;
  for(Real Q:profile?std::vector<Real>{91.1876L}:std::vector<Real>{75.L,91.1876L,120.L})
    for(Real mass_factor:profile?std::vector<Real>{1}:std::vector<Real>{.5L,1.L,2.L})
    for(Real scale_factor:profile?std::vector<Real>{1}:std::vector<Real>{.5L,1.L,2.L})
    for(Real alpha_factor:profile?std::vector<Real>{1}:std::vector<Real>{1.L,.125L,.0625L,.03125L}) {
      QEDParameters q;q.muon_mass*=mass_factor;q.mu_0=q.muon_mass;
      q.mu_frag=Q*scale_factor;q.alpha*=alpha_factor;
      const leptonff_dyturbo::MassRadiator ll(q,grid,PhotonicAccuracy::ll),nll(q,grid);
      const leptonff_dyturbo::MassRadiator quadratic(q,grid,PhotonicAccuracy::nll,
          leptonff_dyturbo::MassMatching::quadratic_finite);
      for(Real N:{1.L,2.L,3.L,5.L,10.L}) {
        const auto l=ll.moment(Q,N),n=nll.moment(Q,N);
        const Real cut=N==1?.25L:N==2?.5L:N==3?.9L:N==5?.99L:.999L;
        const auto bl=ll.bin(Q,Q*std::sqrt(cut),Q*1.01L),bn=nll.bin(Q,Q*std::sqrt(cut),Q*1.01L);
        const auto qm=quadratic.moment(Q,N),qb=quadratic.bin(Q,Q*std::sqrt(cut),Q*1.01L);
        if(!first)std::cout<<',';first=false;
        std::cout<<"{\"Q\":"<<Q<<",\"mass\":"<<q.muon_mass<<",\"scale\":"<<q.mu_frag
          <<",\"alpha\":"<<q.alpha<<",\"N\":"<<N<<",\"order_alpha\":"<<n[1]
          <<",\"ll\":"<<l[2]<<",\"nll\":"<<n[2]<<",\"cut\":"<<cut
          <<",\"bin_order_alpha\":"<<bn[1]<<",\"bin_ll\":"<<bl[2]<<",\"bin_nll\":"<<bn[2]
          <<",\"quadratic\":"<<qm[2]<<",\"bin_quadratic\":"<<qb[2]<<"}";
      }
    }
  std::cout<<"]}\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
