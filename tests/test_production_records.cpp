#include "leptonff/production_records.hpp"
#include "check.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

int main() {
  using namespace leptonff;
  const auto path=std::filesystem::temp_directory_path()/("leptonff-records-"+std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::create_directory(path);
  try {
    ProductionMetadata meta;meta.producer="fixture";meta.physics_id="physical-config";meta.run_id="run-1";meta.random_seed=1045;meta.weight_names={"central","scale"};
    const BornDecayPoint p{{50,40,0,30},{50,-40,0,-30},123};
    ProductionWriter writer((path/"points.lff").string(),meta);
    for(int iteration=1;iteration<=2;++iteration) {
      writer.begin();writer.add({p,{1.2L*iteration,2.4L*iteration}});writer.add({p,{-.3L*iteration,-.6L*iteration}});writer.commit();
      writer.begin();writer.add({p,{.5L*iteration,1.L*iteration}});writer.commit();
      writer.begin();writer.add({p,{1000,2000}});writer.commit(false);
      writer.begin();writer.commit();writer.block({1,1,iteration,4,0});
    }
    writer.finish();
    throws<std::invalid_argument>([&]{ProductionWriter duplicate((path/"points.lff").string(),meta);},"refuse archive overwrite");
    ProductionReader reader((path/"points.lff").string());GroupedIntegration integral(2);ProductionBlock block;
    int count=0;
    while(reader.next_block(block,[&](const ProductionGroup& g){
      std::vector<Real> sum(2);
      for(const auto& t:g.terms) {
        close(t.point.first.px,p.first.px,0,"lossless momentum round trip");
        if(t.point.radiation_seed!=p.radiation_seed)throw std::runtime_error("Seed not preserved");
        for(int j=0;j<2;++j)sum[j]+=t.weights[j];
      }
      integral.add_group(sum);++count;
    }))integral.end_block(block);
    if(count!=4)throw std::runtime_error("Zero/rejected groups were stored");
    const auto result=integral.estimate();
    close(result.mean[0],2.1L,1e-15L,"common-iteration mean, no extra ncall division");
    close(result.covariance[0][0],.95L,1e-14L,"sum real/dipole BEFORE squaring; include zero trials");
    close(result.mean[1],2*result.mean[0],1e-15L,"weight-vector replay");
    close(result.covariance[0][1],2*result.covariance[0][0],1e-14L,"retain correlated variation weights");
    GroupedIntegration resumed(1);
    resumed.add_group({.9L});resumed.add_group({.5L});resumed.end_block({1,1,1,4,2});
    resumed.add_group({.45L});resumed.add_group({.25L});resumed.end_block({1,1,11,8,2});
    const auto continuation=resumed.estimate();
    close(continuation.mean[0],1.05L,1e-15L,"resumed host iterations use common weights");
    close(continuation.covariance[0][0],(.76L+1.63L/7)/4,1e-14L,"resumed blocks retain their own zero-trial denominators");
    std::ifstream input(path/"points.lff");std::ostringstream contents;contents<<input.rdbuf();std::string data=contents.str();
    data.resize(data.size()-4);std::ofstream(path/"interrupted.lff")<<data;
    throws<std::invalid_argument>([&]{ProductionReader r((path/"interrupted.lff").string());while(r.next_block(block,[](const ProductionGroup&){})){};},"reject missing final commit marker");
    auto invalid=meta;invalid.contains_qed=true;
    throws<std::invalid_argument>([&]{ProductionWriter r((path/"bad.lff").string(),invalid);},"reject QED double counting");
    invalid=meta;invalid.full_born_angles=false;
    throws<std::invalid_argument>([&]{ProductionWriter r((path/"bad.lff").string(),invalid);},"reject angular-precut samples");
    throws<std::invalid_argument>([&]{integral.end_block({1,1,1,4,0});},"reject repeated block");
    std::filesystem::remove_all(path);
  }catch(...){std::filesystem::remove_all(path);throw;}
  std::cout<<"Production record round trip, subtraction groups, normalisation and covariance passed\n";
}
