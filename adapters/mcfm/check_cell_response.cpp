// Deterministic grid-response diagnostic, not an integrated cross section.
#include "leptonff/coupled_lab.hpp"
#include <iomanip>
#include <iostream>
#include <cmath>
using namespace leptonff;
int main() {
  std::cout<<"process,bins,measure,point,ff_higher,projected_higher,geometry_higher,total_higher,order2\n"<<std::setprecision(18);
  SpinTensor rho{};for(int i=0;i<3;++i)rho[i][i]=1.L/3;
  for(auto kind:{SpinDecay::wminus,SpinDecay::neutral_vector}) {
    const bool neutral=kind==SpinDecay::neutral_vector;
    const Real M=neutral?91.1876L:80.385L;
    const FourVector parent{std::sqrt(M*M+30*30+50*50),30,0,50};
    const auto born=born_decay(M,.3L,.7L);
    const auto first=rest_frame_boost(born.first,parent,1),second=rest_frame_boost(born.second,parent,1);
    MomentEvolutionParameters p;p.mu_frag=M;
    for(std::size_t bins:{4096,8192,16384,32768})for(auto mode:{FractionMeasure::atoms,FractionMeasure::conservative_cells}) {
      QEDXSpaceGrid grid;grid.bins=bins;
      const CoupledLabDressing d(kind,p,grid,1e-4L,mode);
      if(d.fragmentation().same().measure()!=mode)throw std::runtime_error("Fraction-measure diagnostic mismatch");
      const auto ff=neutral?d.fragmentation().z(first,second):d.fragmentation().w(first,second);
      const auto h=neutral?d.z_projected(first,second):d.w_projected(first,second);
      const auto f=[&](const FourVector& a,const FourVector& b,const FourVector& k) {
        const Real bare=neutral?(FiducialCuts{}.accepts(a,b)?1:0):(a.pt()>25 && std::abs(a.eta())<2.4L?1:0);
        return std::array<Real,3>{bare,neutral?d.z_real(a,b,k):d.w_real(a,b,k),neutral?d.z_real(a,b,k,{},true):d.w_real(a,b,k,{},true)};
      };
      for(auto point:{std::array<Real,2>{.45L,.15L},{.25L,.4L}}) {
        const auto tensor=decay_hard_geometry_channels(kind,M,point[0],point[1],decay_rotation(.45L,.7L,1.1L),parent,f,neutral?0:1);
        const std::array<Real,3> g{spin_contract(rho,tensor[0]),spin_contract(rho,tensor[1]),spin_contract(rho,tensor[2])};
        const auto result=CoupledLabDressing::match(ff,h,g);
        const Real projected=h.a_hard*h.evolved-h.a_initial*h.lattice;
        const Real geometry=h.a_hard*g[1]-h.a_initial*g[0];
        std::cout<<(neutral?"z":"wminus")<<','<<bins<<','<<(mode==FractionMeasure::atoms?"atoms":"cells")<<','<<point[0]<<','
          <<ff.matched-ff.fixed_order<<','<<projected<<','<<geometry<<','<<result.matched-result.fixed_order<<','
          <<h.a_initial*h.a_initial*result.coefficient2<<'\n';
      }
    }
  }
}
