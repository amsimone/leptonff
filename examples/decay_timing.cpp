#include <leptonff/matched.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
using namespace leptonff;
using Clock=std::chrono::steady_clock;
volatile Real checksum=0;
template<class F> double time(unsigned count,const F& f) {
  f();const auto start=Clock::now();for(unsigned i=0;i<count;++i)f();
  return std::chrono::duration<double,std::milli>(Clock::now()-start).count()/count;
}
}
int main(int argc,char** argv) {
  try {
    const auto repetitions=argc>1?std::stoul(argv[1]):256;
    if(!repetitions || repetitions>100000)throw std::invalid_argument("Repetitions must be in [1,100000]");
    std::cout<<"process,model,bins,construct_ms,point_ms,batch4_ms_per_point,response_ms,cached_contract_ms\n";
    for(bool neutral:{false,true})for(auto model:{DecayPreset::photonic_nll,DecayPreset::leptonic_nll}) {
      const Real M=neutral?91.1876L:80.4L;
      const FourVector q{std::sqrt(M*M+1300),20,0,30};
      const BornDecayPoint point{
        rest_frame_boost({M/2,.9L*M/2,0,std::sqrt(.19L)*M/2},q,1),
        rest_frame_boost({M/2,-.9L*M/2,0,-std::sqrt(.19L)*M/2},q,1),42};
      auto c=make_decay_configuration(neutral?SpinDecay::neutral_vector:SpinDecay::wplus,model,
          AlphaInput::on_shell(1/137.035999084L),.1056583755L,M);
      const auto start=Clock::now();const DecayCalculator decay(c);
      const auto setup=std::chrono::duration<double,std::milli>(Clock::now()-start).count();
      for(unsigned subdivision:{1,4}) {
        const MeasurementPlan measurement(neutral,LabSpectrum::configured(neutral,LabSpectrum::reference_binning(neutral,subdivision)));
        const auto point_ms=time(repetitions,[&]{checksum=decay.evaluate(point,measurement).total.matched;});
        std::vector<BornDecayPoint> points(256,point);
        const auto batch_ms=time(2,[&]{checksum=decay.evaluate_batch(points,measurement,4).front().total.matched;})/points.size();
        const auto response_ms=time(2,[&]{checksum=decay.response(q,measurement,42,{8,16}).total[2][0][0].real();});
        const auto cached=decay.response(q,measurement,42,{8,16});
        SpinTensor rho{};rho[0][0]=1;rho[1][1]=.3L;rho[2][2]=-.1L;
        const auto contraction=time(1024,[&]{checksum=cached.contract(rho).total.matched;});
        std::cout<<(neutral?"Z":"W")<<','<<(model==DecayPreset::photonic_nll?"photonic_nll":"leptonic_nll")<<','
          <<measurement.spectrum().bins.size()<<','<<setup<<','<<point_ms<<','<<batch_ms<<','<<response_ms<<','<<contraction<<'\n';
      }
    }
  }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
