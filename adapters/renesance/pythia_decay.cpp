// SPDX-License-Identifier: MIT
#include "pythia_decay.hpp"
#include "Pythia8/Pythia.h"
#include <stdexcept>
namespace leptonff_reference {
struct PythiaDecay::Impl {Pythia8::Pythia p;};
std::string PythiaDecay::version(){return std::to_string(PYTHIA_VERSION);}
PythiaDecay::PythiaDecay(double alpha,double cutoff,bool matrix_element):impl_(new Impl) {
  if(!(alpha>0 && alpha<1 && cutoff>0))throw std::invalid_argument("Invalid PYTHIA inputs");
  auto& p=impl_->p;
  const auto set=[&](const std::string& s){if(!p.readString(s))throw std::runtime_error("PYTHIA rejected "+s);};
  set("Print:quiet = on");set("ProcessLevel:all = off");
  set("PartonLevel:ISR = off");set("PartonLevel:MPI = off");
  set("HadronLevel:all = off");set("TimeShower:QCDshower = off");
  set("TimeShower:QEDshowerByQ = off");set("TimeShower:QEDshowerByL = on");
  set("TimeShower:QEDshowerByGamma = off");
  set("TimeShower:alphaEMorder = 0");
  p.settings.parm("StandardModel:alphaEM0",alpha);
  p.settings.parm("TimeShower:pTminChgL",cutoff);
  if(std::abs(p.settings.parm("StandardModel:alphaEM0")-alpha)>1e-14 ||
     std::abs(p.settings.parm("TimeShower:pTminChgL")-cutoff)>1e-14)
    throw std::runtime_error("PYTHIA clamped the requested coupling or cutoff");
  p.settings.flag("TimeShower:MEcorrections",matrix_element);
  p.settings.flag("TimeShower:MEafterFirst",matrix_element);
  p.particleData.m0(13,.1056583755);
  set("Random:setSeed = on");set("Random:seed = 1");
  if(!p.init())throw std::runtime_error("PYTHIA initialization failed");
}
PythiaDecay::~PythiaDecay()=default;
std::pair<leptonff::FourVector,leptonff::FourVector> PythiaDecay::draw(
    leptonff::SpinDecay kind,const leptonff::BornDecayPoint& point,std::uint64_t seed) {
  using namespace leptonff;
  auto& p=impl_->p;
  const auto a=point.first,b=point.second;
  const FourVector parent{a.energy+b.energy,a.px+b.px,a.py+b.py,a.pz+b.pz};
  const auto M=invariant_mass(a,b);const auto rest=rest_frame_boost(a,parent,-1);
  const Real m=.1056583755L,norm=std::hypot(rest.pt(),rest.pz);
  const bool z=kind==SpinDecay::neutral_vector;
  const int sign=kind==SpinDecay::wplus?1:-1;
  const double E=z?M/2:(M*M+m*m)/(2*M),k=std::sqrt(E*E-m*m);
  p.event.reset();
  const int mother=p.event.append(z?23:24*sign,-22,0,0,0,0,0,0,0,0,0,M,M);
  const int first=p.event.append(z?13:-13*sign,23,mother,0,0,0,0,0,
      k*rest.px/norm,k*rest.py/norm,k*rest.pz/norm,E,m);
  const int second=p.event.append(z?-13:14*sign,23,mother,0,0,0,0,0,
      -k*rest.px/norm,-k*rest.py/norm,-k*rest.pz/norm,M-E,z?m:0);
  p.event[mother].daughters(first,second);
  p.event[first].scale(M/2);p.event[second].scale(M/2);
  p.rndm.init(int(seed%899999999)+1);
  p.forceTimeShower(first,second,M/2);
  const auto read=[&](int i){i=p.event[i].iBotCopyId();const auto& v=p.event[i];
    if(!v.isFinal())throw std::runtime_error("PYTHIA lepton ancestry did not reach the final state");
    return rest_frame_boost({v.e(),v.px(),v.py(),v.pz()},parent,1);};
  FourVector sum{};
  for(int i=1;i<p.event.size();++i)if(p.event[i].isFinal()){
    sum.energy+=p.event[i].e();sum.px+=p.event[i].px();sum.py+=p.event[i].py();sum.pz+=p.event[i].pz();}
  if(std::max({std::abs(sum.energy-M),std::abs(sum.px),std::abs(sum.py),std::abs(sum.pz)})>1e-7L*M)
    throw std::runtime_error("PYTHIA momentum conservation failure");
  return {read(first),read(second)};
}
}
