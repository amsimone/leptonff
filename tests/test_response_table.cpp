#include "leptonff/response_table.hpp"
#include "leptonff/decay_setup.hpp"
#include "check.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

int main() {
  using namespace leptonff;
  const auto path=std::filesystem::temp_directory_path()/("leptonff-response-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(path);
  try {
    for(bool neutral:{false,true}) {
      auto config=make_decay_configuration(neutral?SpinDecay::neutral_vector:SpinDecay::wplus,DecayPreset::photonic_nll,AlphaInput::on_shell(1/137.035999084L));
      config.grid.bins=128;config.hard_samples=2;
      const DecayCalculator c(config);const MeasurementPlan plan(neutral,LabSpectrum::reference(neutral));
      const DecayAngularGrid grid{4,8};
      const std::vector<DecayResponseNode> nodes{{{100,0,0,0},42},{{110,10,0,20},137}};
      const auto table=DecayResponseTable::build(c,plan,nodes,grid,2);
      const auto filename=(path/(neutral?"z.lfr":"w.lfr")).string();table.save(filename);
      const auto loaded=DecayResponseTable::load(filename,c,plan,grid);
      SpinTensor rho{};rho[0][0]=1;rho[1][1]=-.2L;rho[2][2]=.7L;rho[0][1]={.1L,.2L};rho[1][0]=std::conj(rho[0][1]);
      for(std::size_t k=0;k<nodes.size();++k) {
        const auto direct=c.response(nodes[k].parent,plan,nodes[k].radiation_seed,grid).contract(rho);
        const auto replay=loaded.response(k).contract(rho);
        for(int j=0;j<5;++j)close(replay.total.values()[j],direct.total.values()[j],0,"persisted exact-node response vs direct tensor calculation");
        for(std::size_t i=0;i<direct.bins.size();++i)for(int j=0;j<5;++j)close(replay.bins[i].values()[j],direct.bins[i].values()[j],0,"persisted bin coefficients round trip exactly");
      }
      throws<std::invalid_argument>([&]{table.save(filename);},"response overwrite guard");
      throws<std::invalid_argument>([&]{DecayResponseTable::load(filename,c,plan,{8,16});},"angular settings invalidate cache");
      config.soft_cutoff*=2;const DecayCalculator changed(config);
      throws<std::invalid_argument>([&]{DecayResponseTable::load(filename,changed,plan,grid);},"physics/sampling settings invalidate cache");
      const MeasurementPlan finer(neutral,LabSpectrum::configured(neutral,LabSpectrum::reference_binning(neutral,2)));
      throws<std::invalid_argument>([&]{DecayResponseTable::load(filename,c,finer,grid);},"binning invalidates cache");
      std::ifstream input(filename);std::ostringstream data;data<<input.rdbuf();auto bytes=data.str();
      auto end=bytes.rfind("END ");bytes.resize(end);std::ofstream(path/"incomplete")<<bytes;
      throws<std::invalid_argument>([&]{DecayResponseTable::load((path/"incomplete").string(),c,plan,grid);},"incomplete cache rejected");
      bytes=data.str();auto location=bytes.find("NODE 100");bytes.replace(location,8,"NODE 101");std::ofstream(path/"corrupt")<<bytes;
      throws<std::invalid_argument>([&]{DecayResponseTable::load((path/"corrupt").string(),c,plan,grid);},"corrupted finite node rejected by checksum");
    }
    std::filesystem::remove_all(path);
  }catch(...) {std::filesystem::remove_all(path);throw;}
  std::cout<<"Reusable exact-node responses, signed contractions, configuration identity and corruption checks passed\n";
}
