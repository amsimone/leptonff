// Legacy public API benchmark, also compilable against leptonff 0.26.
// Keep its fixed point, seed and controls unchanged for before/after timing.
#include "leptonff/production_records.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace leptonff;
using Clock=std::chrono::steady_clock;
volatile Real sink=0;
template<class F> void bench(const std::string& name,F f,int limit=20000) {
  f();int n=0;auto t=Clock::now();double s;
  do {f();++n;s=std::chrono::duration<double>(Clock::now()-t).count();}while(n<limit && s<.25);
  std::cout<<name<<","<<n<<","<<s*1000/n<<std::endl;
}
int main(int argc,char** argv) {
  std::cout<<std::setprecision(8)<<"# sizeof(Real)="<<sizeof(Real)<<"\nname,repeats,ms_per_call\n";
  const std::string task=argc>1?argv[1]:"events";
  if(task=="stats") {
    for(int n:{90,330,450,1650})for(bool sparse:{false,true}) {
      GroupedIntegration g(n);std::vector<Real> x(n);
      for(int i=0;i<n;++i)x[i]=sparse && i%10?0:Real(i+1)/n;
      bench("cov_"+std::to_string(n)+(sparse?"_10pct":"_dense"),[&]{g.add_group(x);});
    }return 0;
  }
  for(bool neutral:{false,true}) {
    const Real M=neutral?91.1876L:80.4L;
    FourVector q{std::sqrt(M*M+20*20+30*30),20,0,30};
    FourVector a=rest_frame_boost({M/2,.9L*M/2,0,std::sqrt(.19L)*M/2},q,1);
    FourVector b=rest_frame_boost({M/2,-.9L*M/2,0,-std::sqrt(.19L)*M/2},q,1);
    const std::string prefix=neutral?"Z":"W";
    bench(prefix+"_projected_hard",[&]{sink=neutral?z_projected_decay_hard(a,b,M):w_projected_decay_hard(a,b,M);});
    if(task=="hard")continue;
    for(auto mode:{DecayMode::photonic_ll,DecayMode::photonic_nll,DecayMode::coupled_nll}) {
      if(task=="profile" && (!neutral || mode!=DecayMode::photonic_nll))continue;
      if(task=="profile_coupled" && (!neutral || mode!=DecayMode::coupled_nll))continue;
      if(task=="leptonic" && mode!=DecayMode::coupled_nll)continue;
      DecayConfiguration c;c.process=neutral?SpinDecay::neutral_vector:SpinDecay::wplus;c.mode=mode;c.evolution.mu_frag=M;
      if(task=="leptonic") {
        c.evolution.flavour_scheme=QEDFlavourScheme::leptons();
        c.evolution.alpha_msbar_at_mu0=leptonic_alpha_msbar(c.evolution.alpha_msbar_at_mu0,c.evolution.mu_0,*c.evolution.flavour_scheme);
        c.fraction_measure=FractionMeasure::conservative_cells;
      }
      const std::string model=mode==DecayMode::photonic_ll?"LL":mode==DecayMode::photonic_nll?"NLL":"coupled";
      auto start=Clock::now();DecayCalculator calculator(c);
      std::cout<<prefix<<"_"<<model<<"_construct,1,"<<std::chrono::duration<double,std::milli>(Clock::now()-start).count()<<std::endl;
      if(task=="profile" || task=="profile_coupled") {
        auto spectrum=LabSpectrum::configured(neutral,LabSpectrum::reference_binning(neutral,4));
        auto t=Clock::now();std::uint64_t seed=0;
        do {auto r=calculator.evaluate({a,b,++seed},spectrum);sink=r.total.matched;}
        while(std::chrono::duration<double>(Clock::now()-t).count()<12);
        std::cout<<"# profile events="<<seed<<std::endl;return 0;
      }
      for(unsigned sub:{1,4}) {
        auto spectrum=LabSpectrum::configured(neutral,LabSpectrum::reference_binning(neutral,sub));
        bench(prefix+"_"+model+"_bins"+std::to_string(spectrum.bins.size()),[&]{auto r=calculator.evaluate({a,b,42},spectrum);sink=r.total.matched;});
      }
      auto spectrum=LabSpectrum::reference(neutral);
      bench(prefix+"_"+model+"_response_8x16",[&]{auto r=calculator.response(q,spectrum,42,{8,16});sink=r.total[2][0][0].real();},2);
    }
  }
}
