// SPDX-License-Identifier: MIT
// Inclusive neutral-current mass only. No laboratory lepton cuts.
#include "leptonff/photonic_dressing.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc,char** argv) {
  using namespace leptonff;
  try {
    if(argc!=4 && argc!=5)throw std::invalid_argument("Usage: mass_response production.tsv output.tsv ff_bins [alpha]");
    std::ifstream input(argv[1]);if(!input)throw std::runtime_error("Cannot read production");
    if(std::ifstream(argv[2]))throw std::runtime_error("Refusing to overwrite migration output");
    LLGrid grid;grid.bins=std::stoul(argv[3]);grid.log_span=32;
    grid.soft_log_cutoff=grid.spacing()/2;
    QEDParameters q;q.alpha=1/137.035999084L;q.mu_frag=91.1876L;
    if(argc==5)q.alpha=std::stold(argv[4]);
    const PhotonicLabDressing decay(q,grid,1e-5L,PhotonicAccuracy::nll);
    const Real a=q.alpha/(2*std::acos(-1.L));
    // Visible bins [50,51],...,[115,116], plus the complement below 50.
    // The production input is restricted to 66..116 GeV. Photons cannot
    // migrate a pair above its parent mass. Inclusive normalization is one.
    std::vector<std::array<Real,3>> sum(67),variance(67);
    Real M,weight,error;unsigned rows=0;
    while(input>>M>>weight>>error) {
      if(!(M>66 && M<116) || error<0)throw std::invalid_argument("Invalid production node");
      const FourVector p{M/2,M/2,0,0},r{M/2,-M/2,0,0};
      std::array<Real,3> accounted{};
      std::vector<std::array<Real,3>> response(67);
      for(unsigned i=0;i<66;++i) {
        const Real lo=50+i,hi=lo+1;
        // Positive dummy pT support exceeds minimum_z, but cannot restrict
        // any visible mass bin: z1*z2>(50/116)^2 implies each z>0.18.
        const FiducialCuts cuts{M*1e-4L/2,100,lo,hi};
        const auto ff=decay.fragmentation().z(p,r,cuts);
        const Real hard=z_projected_decay_hard(p,r,q.mu_frag,cuts);
        response[i+1]={ff.born,ff.fixed_order+a*hard,
          ff.matched+a*(hard+decay.z_projected_shift(p,r,cuts))};
        for(unsigned j=0;j<3;++j)accounted[j]+=response[i+1][j];
      }
      for(unsigned j=0;j<3;++j)response[0][j]=1-accounted[j];
      for(unsigned i=0;i<67;++i)for(unsigned j=0;j<3;++j) {
        sum[i][j]+=weight*response[i][j];variance[i][j]+=std::pow(error*response[i][j],2);
      }
      ++rows;
    }
    if(!input.eof() || !rows)throw std::runtime_error("Malformed/empty production table");
    std::ofstream out(argv[2]);out<<std::setprecision(17);
    out<<"# lower upper no_qed order_alpha photonic_nll no_qed_error order_alpha_error photonic_nll_error (fb)\n";
    for(unsigned i=0;i<67;++i) {
      out<<(i?49+i:0)<<' '<<(i?50+i:50);
      for(Real x:sum[i])out<<' '<<x;
      for(Real x:variance[i])out<<' '<<std::sqrt(x);
      out<<'\n';
    }
    if(!out)throw std::runtime_error("Cannot write mass spectrum");
    std::cout<<"Convolved "<<rows<<" parent-mass bins with fixed-alpha photonic NLL and finite O(alpha) matching\n";
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
