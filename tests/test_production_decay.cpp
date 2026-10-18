#include "leptonff/production_decay.hpp"
#include "check.hpp"
#include <future>
#include <map>

int main() {
  using namespace leptonff;
  const FourVector mu{50,40,0,30},other{50,-40,0,-30};
  for(auto kind:{SpinDecay::wplus,SpinDecay::wminus,SpinDecay::neutral_vector}) {
    const bool neutral=kind==SpinDecay::neutral_vector;
    DecayConfiguration c;c.process=kind;c.grid.bins=128;c.hard_samples=4;
    c.evolution.mu_frag=neutral?91.1876L:80.385L;
    const DecayCalculator decay(c);
    const auto spectrum=LabSpectrum::reference(neutral);
    const BornDecayPoint point{mu,other,12345};
    const auto p=decay.evaluate(point,spectrum);
    auto job1=std::async(std::launch::async,[&]{return decay.evaluate(point,spectrum);});
    auto job2=std::async(std::launch::async,[&]{return decay.evaluate(point,spectrum);});
    const auto repeated=job1.get();const auto simultaneous=job2.get();
    for(int j=0;j<5;++j)close(p.total.values()[j],repeated.total.values()[j],0,"immutable evaluator repeatability");
    for(int j=0;j<5;++j)close(p.total.values()[j],simultaneous.total.values()[j],0,"concurrent immutable context");
    std::map<std::string,std::array<Real,5>> totals;
    for(std::size_t k=0;k<p.bins.size();++k)for(int j=0;j<5;++j)
      totals[spectrum.bins[k].observable][j]+=p.bins[k].values()[j];
    for(const auto& entry:totals)for(int j=0;j<5;++j)
      close(entry.second[j],p.total.values()[j],2e-12L,"custom-channel spectrum sums to rate");
    close(p.total.matched-p.total.order1,p.total.higher,0,"higher-order difference");
    const auto fine=LabSpectrum::configured(neutral,LabSpectrum::reference_binning(neutral,4));
    if(fine.bins.size()<=3*spectrum.bins.size())throw std::runtime_error("Finer observable bins missing");
    const auto f=decay.evaluate(point,fine);
    close(f.total.matched,p.total.matched,0,"binning-independent rate");
    c.evolution.alpha_msbar_at_mu0=0;const DecayCalculator zero(c);
    const auto tensor=zero.response({100,0,0,0},spectrum,123,{32,16});
    WLabCuts inclusive_w;inclusive_w.muon_pt_min=1e-6L;inclusive_w.muon_eta_max=100;
    FiducialCuts inclusive_z;inclusive_z.pt_min=1e-6L;inclusive_z.eta_max=100;
    const auto inclusive=LabSpectrum::configured(neutral,neutral?std::vector<ObservableBinning>{{"mass",{60,120}}}:
      std::vector<ObservableBinning>{{"pt",{1e-6L,std::numeric_limits<Real>::infinity()}}},inclusive_w,inclusive_z);
    const auto unit=zero.response({100,0,0,0},inclusive,123,{4,8});
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
      close(std::abs(unit.total[0][i][j]-(i==j?1.L:0.L)),0,2e-13L,"inclusive Born tensor integrates to identity");
    SpinTensor rho{};rho[0][0]=1;rho[1][1]=2;rho[2][2]=-.2L;rho[0][1]={.1L,.2L};rho[1][0]=std::conj(rho[0][1]);
    const auto contracted=tensor.contract(rho);
    close(contracted.total.matched,contracted.total.born,1e-13L,"zero-QED tensor response identity");
    close(contracted.total.order1,contracted.total.born,1e-13L,"zero-QED tensor fixed order identity");
    for(int i=0;i<3;++i)rho[i][i]*=-2;
    rho[0][1]*=-2;rho[1][0]*=-2;
    close(tensor.contract(rho).total.born,-2*contracted.total.born,1e-12L,"signed linear production-tensor contraction");
    const auto finite=decay.response({100,0,0,0},spectrum,123,{4,8});
    const auto rate=finite.contract(rho).total;
    // Independent weighted Born-angle route, including the dual projector
    // in each event. Four-point Gauss x eight azimuths integrates B*K exactly.
    const Real nodes[]={-.8611363115940525752L,-.3399810435848562648L,.3399810435848562648L,.8611363115940525752L};
    const Real weights[]={.3478548451374538574L,.6521451548625461426L,.6521451548625461426L,.3478548451374538574L};
    std::array<Real,5> angular{};
    for(int t=0;t<4;++t)for(int k=0;k<8;++k) {
      const Real phi=2*std::acos(-1.L)*(k+.5L)/8,cost=nodes[t],sint=std::sqrt(1-cost*cost);
      const std::array<Real,3> n{sint*std::cos(phi),sint*std::sin(phi),cost};
      const BornDecayPoint node{{50,50*n[0],50*n[1],50*n[2]},{50,-50*n[0],-50*n[1],-50*n[2]},123};
      const auto prediction=decay.cumulative(node,{spectrum.selections[0]})[0].values();
      const Real weight=weights[t]/16*spin_contract(rho,spin_born_tensor(kind,n));
      for(int j=0;j<5;++j)angular[j]+=weight*prediction[j];
    }
    for(int j=0;j<5;++j)close(angular[j],rate.values()[j],2e-11L,"direct tensor response vs independent angular tomography");
    close(rate.matched-rate.order1,rate.higher,1e-12L,"finite response remains linear after matching");
    std::map<std::string,Real> sum;
    const auto bins=finite.contract(rho);
    for(std::size_t k=0;k<bins.bins.size();++k)sum[spectrum.bins[k].observable]+=bins.bins[k].matched;
    for(auto entry:sum)close(entry.second,rate.matched,2e-12L,"tensor response observable coverage");
  }
  for(auto kind:{SpinDecay::wplus,SpinDecay::wminus,SpinDecay::neutral_vector}) {
    DecayConfiguration c;c.process=kind;c.grid.bins=128;c.hard_samples=2;c.mode=DecayMode::coupled_nll;
    const DecayCalculator coupled(c);
    const auto s=LabSpectrum::reference(kind==SpinDecay::neutral_vector);
    const auto value=coupled.response({100,0,0,0},s,18,{4,8});
    if(value.bins.size()!=17)throw std::runtime_error("Coupled tensor response failed");
  }
  const double raw[4]={40,0,30,50};const auto p=map_massless_parent(raw);
  close(p.energy,50,0,"explicit host momentum convention");
  const auto seed=radiation_seed(raw,raw);
  if(seed!=radiation_seed(raw,raw))throw std::runtime_error("Seed instability");
  throws<std::invalid_argument>([]{DecayConfiguration c;c.hard_samples=0;DecayCalculator d(c);},"reject unmatched NLL");
  throws<std::invalid_argument>([]{LabSpectrum::configured(false,{{"pt",{25,30,60}}});},"require overflow coverage");
  throws<std::invalid_argument>([]{LabSpectrum::configured(false,{{"mass",{60,120}}});},"reject wrong-process observable");
  std::cout<<"Generator-independent decay, tensor response, signed weights and refined bins passed\n";
}
