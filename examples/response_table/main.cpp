#include <leptonff/matched.hpp>
#include <iostream>

int main(int argc,char** argv) {
  using namespace leptonff;
  try {
    auto settings=make_decay_configuration(SpinDecay::wplus,DecayPreset::photonic_nll,
        AlphaInput::on_shell(1/137.035999084L),.1056583755L,80.385L);
    settings.grid.bins=1024;settings.hard_samples=16;
    const DecayCalculator decay(settings);
    const MeasurementPlan measurement(false,LabSpectrum::reference(false));
    const std::vector<DecayResponseNode> nodes{{{100,0,0,0},42},{{110,10,0,20},137}};
    const DecayAngularGrid angular{16,32};
    auto table=DecayResponseTable::build(decay,measurement,nodes,angular,2);
    if(argc==2) {table.save(argv[1]);table=DecayResponseTable::load(argv[1],decay,measurement,angular);}
    SpinTensor rho{};rho[0][0]=1;rho[1][1]=.3L;rho[2][2]=-.1L;
    for(std::size_t i=0;i<table.size();++i)
      std::cout<<"node "<<i<<", matched contraction "<<double(table.response(i).contract(rho).total.matched)<<'\n';
    std::cout<<"Illustrative signed production tensor; no cross-section or quadrature-precision claim.\n";
  }catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
