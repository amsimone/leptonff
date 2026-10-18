// SPDX-License-Identifier: MIT
// Standalone fixed-Q muon FSR benchmark; upstream KKMC supplies all amplitudes.
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
extern "C" {
void kk2f_readatax_(const char*,int*,int*,double*,std::size_t);
void kk2f_initialize_(double*);
void kk2f_make_();
void kk2f_getwtlist_(double*,double*);
void kk2f_getfermions_(double*,double*);
void kk2f_getphotall_(int*,double*);
void pseumar_initialize_(int*,int*,int*);
}
namespace {
constexpr unsigned modes=5,dim=2*modes;
struct Observable {
  std::string name;double value;
  std::array<long double,dim> sum{};
  std::array<long double,dim*dim> product{};
  void add(const std::array<double,modes>& w,double measurement) {
    std::array<long double,dim> x;
    for(unsigned k=0;k<modes;++k){x[k]=w[k];x[k+modes]=w[k]*measurement;}
    for(unsigned k=0;k<dim;++k){sum[k]+=x[k];for(unsigned j=0;j<dim;++j)product[k*dim+j]+=x[k]*x[j];}
  }
};
}
int main(int argc,char** argv){try {
  if(argc!=10)throw std::runtime_error("Usage: reference defaults output Q mass alpha trials seed cutoff fsr(0|1)");
  const std::string defaults=argv[1],output=argv[2];
  if(std::ifstream(output))throw std::runtime_error("Refusing to overwrite result");
  const double Q=std::stod(argv[3]),mass=std::stod(argv[4]),alpha=std::stod(argv[5]),cutoff=std::stod(argv[8]);
  const auto trials=std::stoull(argv[6]);int seed=std::stoi(argv[7]),fsr=std::stoi(argv[9]);
  if(!(Q>10*mass && mass>0 && alpha>0 && alpha<.1 && cutoff>0 && cutoff<.001) ||
      trials<2 || seed<=0 || seed>=900000000 || (fsr!=0 && fsr!=1))
    throw std::runtime_error("Invalid fixed-Q benchmark inputs");
  int reset=1,size=10000,zero=0;double x[10000]{};
  kk2f_readatax_(defaults.c_str(),&reset,&size,x,defaults.size());
  const auto set=[&](unsigned k,double v){x[k-1]=v;};
  set(1,Q);set(2,0);set(5,0);set(6,0);set(7,0);
  set(10,1);set(12,0);set(13,0);set(16,cutoff);set(17,1);
  set(20,0);set(21,fsr);set(22,1);set(27,0);set(28,fsr?1:0);set(30,1/alpha);
  set(50,0);set(53,0);set(100,0);
  for(unsigned k=401;k<=420;++k)set(k,0);
  set(413,1);set(636,mass);set(2001,-1);set(2002,-1);
  pseumar_initialize_(&seed,&zero,&zero);
  kk2f_initialize_(x);
  std::vector<Observable> obs;
  for(double N:{2.,3.,5.,10.})obs.push_back({"moment",N});
  for(double cut:{.25,.5,.75,.9,.95,.98,.99,.995,.999})obs.push_back({"cumulative",cut});
  const auto start=std::chrono::steady_clock::now();
  double closure=0;unsigned long long bad_kinematics=0;
  for(unsigned long long event=0;event<trials;++event) {
    kk2f_make_();
    double main_weight=0,weights[1000]{},a[4],b[4],photons[400]{};
    kk2f_getwtlist_(&main_weight,weights);kk2f_getfermions_(a,b);
    int count=0;kk2f_getphotall_(&count,photons);
    std::array<double,modes> w;
    for(unsigned k=0;k<modes;++k)w[k]=weights[std::array<unsigned,modes>{72,73,251,252,253}[k]-1];
    for(auto v:w)if(!std::isfinite(v))throw std::runtime_error("Nonfinite KKMC order weight");
    bool empty=true;for(auto v:w)empty=empty && v==0;
    if(empty)continue; // A null weighted trial need not carry physical momenta.
    double square=0,delta=0,missing_energy=0,missing_space2=0;
    for(unsigned k=0;k<4;++k) {
      const double p=a[k]+b[k];square+=(k==3?1:-1)*p*p;
      double total=p;
      for(int j=0;j<count;++j)total+=photons[k*100+j];
      delta=std::max(delta,std::abs(total-(k==3?Q:0))/Q);
      if(k==3)missing_energy=Q-total;else missing_space2+=total*total;
    }
    closure=std::max(closure,delta);
    const double z=square/(Q*Q);
    // CEEX requires KeyPia=1: photons below Q*cutoff/2 are integrated out
    // of the event record, while their recoil remains in the muon momenta.
    // Bound the missing soft four-vector; do not call this exact closure.
    if(z<0 || z>1+1e-9 || delta>50*cutoff+1e-9 ||
       missing_energy < -1e-9*Q || std::sqrt(missing_space2)>missing_energy+1e-9*Q)
      ++bad_kinematics;
    if(bad_kinematics){
      std::cerr<<std::setprecision(18)<<"event="<<event<<", z="<<z<<", closure="<<delta<<", photons="<<count<<'\n';
      for(unsigned k=0;k<4;++k){std::cerr<<a[k]<<' '<<b[k]<<' ';
        for(int j=0;j<count;++j)std::cerr<<photons[k*100+j]<<' ';std::cerr<<'\n';}
      throw std::runtime_error("Invalid KKMC momenta");
    }
    for(auto& o:obs)o.add(w,o.name=="moment"?std::pow(z,o.value-1):(z>=o.value?1:0));
    if((event+1)%100000==0)std::cerr<<"Completed "<<event+1<<"/"<<trials<<'\n';
  }
  std::ofstream out(output);out<<std::setprecision(18);
  out<<"{\"schema\":\"leptonff-kkmcee-fixed-q-v1\",\"Q\":"<<Q<<",\"mass\":"<<mass
    <<",\"alpha\":"<<alpha<<",\"trials\":"<<trials<<",\"seed\":"<<seed<<",\"cutoff\":"<<cutoff
    <<",\"isr\":false,\"ifi\":false,\"weak\":false,\"pairs\":false,\"running_alpha\":false,\"fsr\":"<<(fsr?"true":"false")
    <<",\"ceex_active\":"<<(fsr?"true":"false")
    <<",\"ceex_min_mass\":"<<x[637]
    <<",\"remove_soft_photons\":true"
    <<",\"normalization\":\"each perturbative order divided by its own inclusive weight\""
    <<",\"modes\":[\"eex1\",\"eex2\",\"ceex0\",\"ceex1\",\"ceex2\"],\"weight_indices\":[72,73,251,252,253]"
    <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
    <<",\"max_resolved_record_missing_fraction\":"<<closure
    <<",\"resolved_record_missing_bound\":"<<50*cutoff+1e-9<<",\"observables\":[";
  bool first=true;
  for(const auto& o:obs){
    if(!first)out<<',';first=false;
    out<<"{\"observable\":\""<<o.name<<"\",\"value\":"<<o.value<<",\"mean\":[";
    for(unsigned k=0;k<dim;++k){if(k)out<<',';out<<o.sum[k]/trials;}
    out<<"],\"covariance\":[";
    for(unsigned k=0;k<dim;++k)for(unsigned j=0;j<dim;++j){
      if(k||j)out<<',';
      out<<(o.product[k*dim+j]-o.sum[k]*o.sum[j]/trials)/(trials*(trials-1.L));
    }
    out<<"]}";
  }
  out<<"]}\n";if(!out)throw std::runtime_error("Failed writing KKMC result");
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
