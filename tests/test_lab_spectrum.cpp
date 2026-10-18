#include "leptonff/lab_spectrum.hpp"
#include "leptonff/coupled_lab.hpp"
#include "check.hpp"
#include <map>

int main() {
  using namespace leptonff;
  QEDXSpaceGrid grid;grid.bins=256;
  for(bool neutral:{false,true}) {
    const auto spectrum=LabSpectrum::reference(neutral);
    if(spectrum.bins.size()!=17)throw std::runtime_error("Incomplete spectrum");
    const FourVector a{65,39,0,52},b{50,-30,40,0};
    std::vector<Real> cumulative;
    for(const auto& s:spectrum.selections)cumulative.push_back(s.accepts(neutral,a,b));
    std::map<std::string,Real> totals;
    const auto value=[&](const std::string& obs,const FourVector& u,const FourVector& v) {
      if(obs=="mass")return invariant_mass(u,v);
      if(obs=="minpt")return std::min(u.pt(),v.pt());
      if(obs=="maxeta")return std::max(std::abs(u.eta()),std::abs(v.eta()));
      if(obs=="pt")return u.pt();
      if(obs=="mt")return transverse_mass(u,v);
      return std::abs(u.eta());
    };
    for(std::size_t i=0;i<spectrum.bins.size();++i) {
      const auto& bin=spectrum.bins[i];const Real x=value(bin.observable,a,b);
      close(spectrum.bin(i,cumulative),cumulative[0]*(x>bin.lower && x<bin.upper?1:0),0,"observed bin vs direct kinematics");
      totals[bin.observable]+=spectrum.bin(i,cumulative);
    }
    for(auto [key,total]:totals)close(total,cumulative[0],0,"all observable bins include overflow");
    // Independent literal signed lattice action on a recoil configuration.
    MomentEvolutionParameters p;p.mu_frag=neutral?91.1876L:80.385L;
    const QEDXSpaceEvolution evolution(p,grid);
    const LabFragmentation ff(evolution);
    for(std::size_t k=0;k<spectrum.bins.size();++k) {
      const auto& bin=spectrum.bins[k];Real literal=0;
      for(std::size_t i=0;i<grid.bins && evolution.z(i)*a.pt()>25;++i)
        for(std::size_t j=0;j<(neutral?grid.bins:1) && evolution.z(j)*b.pt()>(neutral?25:0);++j) {
          const auto u=a.scaled(evolution.z(i)),v=neutral?b.scaled(evolution.z(j)):b;
          const Real x=value(bin.observable,u,v);
          if(!spectrum.selections[0].accepts(neutral,u,v) || x<=bin.lower || x>=bin.upper)continue;
          const auto& w=evolution.weights();
          literal+=neutral?(w[0][0][i]*w[0][0][j]+w[0][1][i]*w[0][1][j])/(evolution.z(i)*evolution.z(j)):w[0][0][i]/evolution.z(i);
        }
      std::vector<Real> evolved;
      for(const auto& s:spectrum.selections)evolved.push_back(neutral?ff.z(a,b,s.z).resummed:ff.w(a,b,s.w).resummed);
      close(spectrum.bin(k,evolved),literal,2e-13L,"differential FF vs direct signed tensor bin");
    }
    // A nontrivial linear measurement and its duplicate share all real and
    // sector currents. Check against the original single-measurement API.
    const auto kind=neutral?SpinDecay::neutral_vector:SpinDecay::wminus;
    const ParentMeasurements f=[](const FourVector& a,const FourVector& b,const FourVector& g) {
      return std::array<Real,3>{a.pt()/100,(b.pt()+g.pt())/100,a.energy/100};
    };
    const FourVector parent{100,0,0,0};const auto rotation=decay_rotation(.3L,.7L,1.2L);
    const auto batch=decay_hard_geometry_channels_batch(kind,100,.7L,.2L,rotation,parent,{f,f},{.4L,.4L});
    const DecayMeasurement bare=[](const FourVector& a,const FourVector& b){return (a.pt()+b.pt())>60?1.L:0.L;};
    const DecayMeasurement evolved=[](const FourVector& a,const FourVector& b){return (a.pt()*a.pt()+b.pt()*b.pt())/10000;};
    const ParentMeasurements photonic=[&](const FourVector& a,const FourVector& b,const FourVector&) {
      return std::array<Real,3>{bare(a,b),evolved(a,b)-bare(a,b),0};
    };
    const auto old_pair=decay_hard_geometry_pair(kind,100,.7L,.2L,rotation,parent,bare,evolved,.4L);
    const auto linear=decay_hard_geometry_linear_batch(kind,100,.7L,.2L,rotation,parent,{photonic},{{.4L,0,0}});
    for(int k=0;k<2;++k)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
      close(std::abs(linear[0][k][i][j]-old_pair[k][i][j]),0,1e-14L,"photonic batch vs independent original difference operator");
    const auto single=decay_hard_geometry_channels(kind,100,.7L,.2L,rotation,parent,f,.4L);
    for(int k=0;k<3;++k)for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
      close(std::abs(batch[0][k][i][j]-single[k][i][j]),0,0,"batch and single geometry identity");
      close(std::abs(batch[0][k][i][j]-batch[1][k][i][j]),0,0,"identical measurements share geometry");
    }
  }
  close(transverse_mass({50,30,0,0},{40,-20,0,0}),std::sqrt(2400.L),1e-14L,"pT-based transverse mass accepts massive observed muons");
  throws<std::invalid_argument>([]{auto c=WLabCuts{};c.muon_pt_min=30;LabSpectrum::reference(false,c);},"changed fiducial edges rejected");
  std::cout<<"Physical differential bins, independent tensor actions and shared geometry passed\n";
}
