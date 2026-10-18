#include <leptonff/matched.hpp>
#include <iostream>

int main() {
  using namespace leptonff;
  DecayConfiguration config;config.process=SpinDecay::wplus;
  config.evolution.mu_frag=80.385L;config.grid.bins=128;config.hard_samples=4;
  const DecayCalculator decay(config);
  const auto spectrum=LabSpectrum::reference(false);
  const auto value=decay.evaluate({{50,40,0,30},{50,-40,0,-30},42},spectrum);
  std::cout<<"Born acceptance "<<value.total.born<<"; matched response "<<value.total.matched<<'\n';
  // A host supplies UNNORMALISED angular coefficients for each production
  // node/matching term. Check its angular/recoil convention before mapping.
  const FourVector parent{125,30,20,60};
  DrellYanAngularWeights production;production.unpolarised=1;production.angular[0]=.2L;
  const auto rho=production_spin_tensor(config.process,production,collins_soper_rotation(parent));
  const auto response=decay.response(parent,spectrum,42,{4,8});
  std::cout<<"Production-tensor contraction "<<response.contract(rho).total.matched<<'\n';
}
