// SPDX-License-Identifier: MIT
// Differential NLO+NLL QCD, with the same linear decay functional in RES/CT/VJ.
#include "dyturbo.h"
#include "settings.h"
#include "resintegr.h"
#include "ctintegr.h"
#include "finintegr.h"
#include "resconst.h"
#include "phasespace.h"
#include "cubacall.h"
#include "leptonff/decay_setup.hpp"
#include "leptonff/production_records.hpp"
#include "photos_decay.hpp"
#include <cuba.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cstring>
#include <chrono>

namespace {
using namespace leptonff;
std::unique_ptr<DecayCalculator> calculator;
std::unique_ptr<MeasurementPlan> measurement;
std::unique_ptr<leptonff_dyturbo::PhotosDecay> photos;
std::unique_ptr<ProductionWriter> records;
int record_stage=-1;
std::uint64_t record_groups=0;
std::array<std::vector<Real>,4> sums,products;
std::array<std::uint64_t,4> trials{};
std::array<Real,4> adaptation_sum{},adaptation_square{};
std::uint64_t retained_trials=0,nonzero=0;
int component=0,photons=1;
bool evaluate_decay=true;
std::string failure;
using Callback=integrand_t (*)(const int&,const double*,const int&,double*,void*,const int&,const int&,double&,const int&);
Callback qcd;
unsigned long long hash_point(const double* x,int n,int iteration) {
  unsigned long long h=1469598103934665603ULL ^ unsigned(opts.rseed) ^ (unsigned(component)<<24) ^ unsigned(iteration);
  for(int k=0;k<n;++k){std::uint64_t b;std::memcpy(&b,x+k,8);h=(h^b)*1099511628211ULL;}
  return h;
}
int integrand(const int& ndim,const double* x,const int& nc,double* f,void* userdata,
    const int& nvec,const int& core,double& weight,const int& iteration) {
  try {
    if(ndim!=(component==0?6:component==1?8:7))throw std::runtime_error("Wrong native integrand dimension");
    if(qcd(ndim,x,nc,f,userdata,nvec,core,weight,iteration)!=nullptr)throw std::runtime_error("Native integrand failure");
    if(!std::isfinite(f[0]) || !std::isfinite(weight))throw std::runtime_error("Non-finite native QCD weight");
    if(iteration<1 || iteration>4)throw std::runtime_error("Unexpected VEGAS iteration");
    const auto stage=std::size_t(iteration-1);
    if(records && record_stage!=int(stage)) {
      if(record_stage>=0)records->block({component,record_stage,1,trials[record_stage],record_groups});
      record_stage=int(stage);record_groups=0;
    }
    ++retained_trials;++trials[stage];if(f[0]==0)return 0;++nonzero;
    const auto& cfg=calculator->configuration();
    const auto first=map_massless_parent(opts.nproc==1?phasespace::p4:phasespace::p3);
    const auto second=map_massless_parent(opts.nproc==1?phasespace::p3:phasespace::p4);
    const BornDecayPoint point{first,second,hash_point(x,ndim,iteration)};
    if(records) {
      // Native weights include 1/N_iteration. Treat each iteration as an
      // additive stratum with the deterministic trial-count fraction.
      records->begin();records->add({point,{Real(weight)*f[0]*iteration/10}});
      records->commit();++record_groups;
    }
    const auto& spectrum=measurement->spectrum();
    std::vector<std::array<Real,5>> values(spectrum.bins.size()+1);
    if(evaluate_decay) {
      const auto prediction=calculator->evaluate(point,*measurement);
      for(std::size_t i=0;i<values.size();++i){const auto& p=i?prediction.bins[i-1]:prediction.total;values[i]={p.born,p.order1,p.matched,0,0};}
    } else {
      std::vector<Real> c;for(const auto& s:spectrum.selections)c.push_back(s.accepts(false,first,second));
      for(std::size_t i=0;i<values.size();++i){const auto b=i?spectrum.bin(i-1,c):c[0];values[i]={b,b,b,0,0};}
    }
    for(int k=0;k<photons;++k) {
      const auto pair=photos->draw(cfg.process,point,point.radiation_seed+0x9e3779b97f4a7c15ULL*unsigned(k+1));
      std::vector<Real> passed;for(const auto& s:spectrum.selections)passed.push_back(s.accepts(false,pair.first,pair.second));
      // PHOTOS-off massive map has the same direction; make it explicitly.
      const FourVector parent{first.energy+second.energy,first.px+second.px,first.py+second.py,first.pz+second.pz};
      const auto rest=rest_frame_boost(first,parent,-1);const auto M=invariant_mass(first,second),m=cfg.evolution.muon_mass;
      const auto E=(M*M+m*m)/(2*M),p=(M*M-m*m)/(2*M),norm=std::hypot(rest.pt(),rest.pz);
      const auto born=rest_frame_boost({E,p*rest.px/norm,p*rest.py/norm,p*rest.pz/norm},parent,1);
      const FourVector nu{parent.energy-born.energy,parent.px-born.px,parent.py-born.py,parent.pz-born.pz};
      std::vector<Real> off;for(const auto& s:spectrum.selections)off.push_back(s.accepts(false,born,nu));
      for(std::size_t i=0;i<values.size();++i){values[i][3]+=(i?spectrum.bin(i-1,off):off[0])/photons;values[i][4]+=(i?spectrum.bin(i-1,passed):passed[0])/photons;}
    }
    const Real w=Real(weight)*f[0];
    for(std::size_t bin=0;bin<values.size();++bin)for(int a=0;a<5;++a) {
      const Real v=w*values[bin][a];sums[stage][5*bin+a]+=v;
      for(int b=0;b<5;++b)products[stage][25*bin+5*a+b]+=v*w*values[bin][b];
    }
    // Train on the fiducial region with a positive floor: finite angular
    // matching also receives contributions from Born-forbidden directions.
    f[0]*=.2+.8*values[0][0];
    const Real adapted=Real(weight)*f[0];
    adaptation_sum[stage]+=adapted;adaptation_square[stage]+=adapted*adapted;
    return 0;
  }catch(const std::exception& e){failure=e.what();return -999;}
}
}
int main(int argc,char** argv) {
  try {
    const char* out=std::getenv("LEPTONFF_DYTURBO_RESULT");if(!out || std::ifstream(out))throw std::runtime_error("Set a new result filename");
    const auto start=std::chrono::steady_clock::now();
    // Use the public initialization stages to reject unsupported scalar
    // callbacks before expensive backend preparation (or extra components).
    DYTurbo::init_const();opts.parse_options(argc,argv);
    if(opts.helicity!=-1 || !opts.resumcpp || !opts.ctcpp || opts.makehistos)
      throw std::runtime_error("Requires unpolarized native C++ routines and client-owned histograms");
    if(!opts.doBORN || !opts.doCT || !opts.doVJ)
      throw std::runtime_error("Requires all RES, CT and V+jet contributions");
    DYTurbo::init_params();
    if(opts.fixedorder || opts.order!=1 || (opts.nproc!=1 && opts.nproc!=2) || opts.makecuts || opts.qed || opts.mellin1d || opts.totpdf!=1)
      throw std::runtime_error("Requires differential NLO+NLL W, no native lepton cuts/QED, one PDF");
    if(opts.order_ct!=1 || opts.order_vjet!=1 || opts.order_hcoef!=1 || opts.order_sudak!=1 || opts.order_evol!=1 || opts.order_expc!=1)
      throw std::runtime_error("Incorrect QCD matching order");
    if(opts.intDimRes!=4 || opts.intDimCT!=6 || opts.intDimVJ!=7 || opts.finord ||
       opts.ptbinwidth || opts.ybinwidth || opts.mbinwidth || opts.dsdqt2 || opts.edsdp3)
      throw std::runtime_error("Requires separate differential RES/CT/VJ with integrated weights");
    const auto setting=[](const char* key,double fallback){const char* x=std::getenv(key);return x?std::stod(x):fallback;};
    auto cfg=make_decay_configuration(opts.nproc==1?SpinDecay::wplus:SpinDecay::wminus,DecayPreset::photonic_nll,
      AlphaInput::on_shell(setting("LEPTONFF_ALPHA",1/137.035999084)),.1056583755L,80.385L);
    cfg.grid.bins=setting("LEPTONFF_GRID_BINS",8192);cfg.hard_samples=setting("LEPTONFF_HARD_SAMPLES",8);
    evaluate_decay=setting("LEPTONFF_EVALUATE_DECAY",1)!=0;photons=int(setting("LEPTONFF_PHOTOS_SAMPLES",1));
    if(photons<1 || photons>1024)throw std::runtime_error("Invalid photon sample count");
    calculator=std::make_unique<DecayCalculator>(cfg);
    std::vector<Real> edges;for(int k=0;k<=48;++k)edges.push_back(2.4L*k/48);
    measurement=std::make_unique<MeasurementPlan>(false,LabSpectrum::configured(false,{{"abseta",edges}}));
    photos=std::make_unique<leptonff_dyturbo::PhotosDecay>(cfg.evolution.alpha_msbar_at_mu0,setting("LEPTONFF_PHOTOS_CUTOFF",1e-5),setting("LEPTONFF_PHOTOS_EXP_EPS",1e-4));
    const auto bins=measurement->spectrum().bins.size()+1;
    std::vector<Real> total(5*bins),covariance(25*bins);
    DYTurbo::WarmUp();DYTurbo::BoundIterator bounds;
    if(bounds.IsEnd())throw std::runtime_error("Empty production domain");
    auto next=bounds;++next;if(!next.IsEnd())throw std::runtime_error("Differential client requires one production-domain bin");
    DYTurbo::SetBounds(bounds);
    if(phasespace::lpt2min!=0 || phasespace::lpt2max<opts.sroot*opts.sroot ||
       phasespace::lexp2ymin>std::exp(-200.) || phasespace::lexp2ymax<std::exp(200.))
      throw std::runtime_error("Native lepton binning must retain full Born support");
    if(!DYTurbo::HasOnlyVegas)throw std::runtime_error("Requires VEGAS for every QCD contribution");
    std::ostringstream details;details<<std::setprecision(17)<<'[';
    const int calls=int(setting("LEPTONFF_CALLS",10000));if(calls<100 || calls%10)throw std::runtime_error("Calls must be >=100 and divisible by 10");
    if(const char* path=std::getenv("LEPTONFF_RECORDS")) {
      const char* identity=std::getenv("LEPTONFF_PRODUCTION_ID");
      if(!identity || !*identity)throw std::runtime_error("Archive export requires LEPTONFF_PRODUCTION_ID with run provenance");
      ProductionMetadata metadata;metadata.process=cfg.process;metadata.producer="DYTurbo 1.4.2 NLO+NLL differential RES/CT/VJ";
      metadata.physics_id=identity;metadata.run_id=out;metadata.random_seed=opts.rseed;
      records=std::make_unique<ProductionWriter>(path,metadata);
    }
    for(component=0;component<3;++component) {
      if(component==0){if(opts.cmw)resconst::initcmw();else resconst::init();qcd=resintegrandMC;}
      else {if(opts.cmwfo)resconst::initcmw();else resconst::init();qcd=component==1?ctintegrandMC:vjlointegrandMC;}
      // intDimCT=6 selects an eight-dimensional C++ MC implementation in
      // 1.4.2: its last coordinates are the two PDF convolutions.
      const int ndim=component==0?6:component==1?8:7;
      for(int stage=0;stage<4;++stage){sums[stage].assign(5*bins,0);products[stage].assign(25*bins,0);}
      trials.fill(0);adaptation_sum.fill(0);adaptation_square.fill(0);retained_trials=nonzero=0;
      int neval=0,fail=0;double value=0,error=0,prob=0;
      Vegas(ndim,1,(integrand_t)integrand,nullptr,1,0.,0.,4,opts.rseed+13*component,calls,calls,
        calls/10,calls/10,1000,0,"",nullptr,&neval,&fail,&value,&error,&prob);
      if(!failure.empty())throw std::runtime_error(failure);
      if(retained_trials!=std::uint64_t(calls))throw std::runtime_error("Unexpected total VEGAS sample count");
      for(int stage=0;stage<4;++stage)if(trials[stage]!=std::uint64_t(calls/10*(stage+1)))throw std::runtime_error("Unexpected iteration sample count");
      if(records){records->block({component,record_stage,1,trials[record_stage],record_groups});record_stage=-1;record_groups=0;}
      if(std::abs(adaptation_sum[3]-value)>1e-9*std::max(1.,std::abs(value)))throw std::runtime_error("Cuba weight normalization mismatch");
      const Real last_var=Real(trials[3])/(trials[3]-1)*(adaptation_square[3]-adaptation_sum[3]*adaptation_sum[3]/trials[3]);
      if(std::abs(last_var-error*error)>1e-7*std::max(1.,error*error))throw std::runtime_error("Cuba last-iteration covariance mismatch");
      if(component)details<<',';
      details<<"{\"part\":"<<component<<",\"evaluations\":"<<neval<<",\"retained\":"<<retained_trials
        <<",\"nonzero\":"<<nonzero<<",\"fixed_budget_fail_flag\":"<<fail<<",\"adaptation_integral_fb\":"<<value<<",\"error_fb\":"<<error<<'}';
      std::cout<<"QCD part "<<component<<": "<<value<<" +/- "<<error<<" fb, retained "<<retained_trials<<std::endl;
      // Deterministic trial-count weights preserve unbiasedness across
      // adaptive iterations. Conditional covariances add as martingale
      // increments; no inverse-variance weights estimated from the outputs.
      for(int stage=0;stage<4;++stage) {
        const Real n=trials[stage],fraction=n/retained_trials;
        for(std::size_t i=0;i<total.size();++i)total[i]+=fraction*sums[stage][i];
        for(std::size_t bin=0;bin<bins;++bin)for(int a=0;a<5;++a)for(int b=0;b<5;++b)
          covariance[25*bin+5*a+b]+=fraction*fraction*n/(n-1)*(products[stage][25*bin+5*a+b]-sums[stage][5*bin+a]*sums[stage][5*bin+b]/n);
      }
    }
    details<<']';if(records)records->finish();DYTurbo::Terminate();
    for(Real value:total)if(!std::isfinite(value))throw std::runtime_error("Non-finite prediction");
    for(Real value:covariance)if(!std::isfinite(value))throw std::runtime_error("Non-finite covariance");
    std::ofstream file(out);file<<std::setprecision(18)<<"{\"schema\":\"leptonff-dyturbo-decay-1\",\"scope\":\"NLO+NLL QCD with pointwise decay in RES, CT and VJ; per-bin covariance\",\"process\":"<<opts.nproc
      <<",\"units\":\"fb\",\"seed\":"<<opts.rseed<<",\"grid\":"<<cfg.grid.bins<<",\"hard_samples\":"<<cfg.hard_samples
      <<",\"observable\":\"abseta\",\"labels\":[\"qed_off\",\"order_alpha\",\"photonic_nll\",\"photos_off\",\"photos_on\"]"
      <<",\"cuts\":{\"muon_pt_min\":25,\"muon_abs_eta_max\":2.4,\"neutrino_pt_min\":0,\"transverse_mass_min\":0}"
      <<",\"qcd\":{\"backend\":\"DYTurbo 1.4.2\",\"order\":\"NLO+NLL with hard coefficient at order one\",\"pdf\":"<<std::quoted(opts.LHAPDFset)
      <<",\"pdf_member\":"<<opts.LHAPDFmember<<",\"sqrts\":"<<opts.sroot<<",\"mass\":["<<phasespace::mmin<<','<<phasespace::mmax
      <<"],\"qt\":["<<phasespace::qtmin<<','<<phasespace::qtmax<<"],\"rapidity\":["<<phasespace::ymin<<','<<phasespace::ymax<<"]}"
      <<",\"photos_samples\":"<<photons<<",\"photos_version\":\""<<leptonff_dyturbo::PhotosDecay::version()<<"\",\"photos_full_me\":false"
      <<",\"photos_cutoff\":"<<setting("LEPTONFF_PHOTOS_CUTOFF",1e-5)<<",\"photos_exponentiation_precision\":"<<setting("LEPTONFF_PHOTOS_EXP_EPS",1e-4)
      <<",\"decay_configuration\":"<<decay_configuration_json(cfg)
      <<",\"integration\":\"four adaptive iterations, deterministic trial-count weighting\",\"parts\":"<<details.str()<<",\"bins\":[";
    for(std::size_t bin=0;bin<bins;++bin){if(bin)file<<',';file<<"{\"lower\":"<<(bin?(bin-1)*.05:0)<<",\"upper\":"<<(bin?bin*.05:2.4)<<",\"mean\":[";
      for(int a=0;a<5;++a)file<<(a?",":"")<<total[5*bin+a];file<<"],\"covariance\":[";
      for(int k=0;k<25;++k)file<<(k?",":"")<<covariance[25*bin+k];file<<"]}";}
    file<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
  }catch(const std::exception& e){std::cerr<<"DYTurbo decay: "<<e.what()<<'\n';return 2;}
}
