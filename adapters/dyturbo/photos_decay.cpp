// SPDX-License-Identifier: MIT
#include "photos_decay.hpp"
#include "Photos.h"
#include "PhotosHEPEVTEvent.h"
#include "PhotosHEPEVTParticle.h"
#include "PH_HEPEVT_Interface.h"
#include <random>
#include <stdexcept>
namespace leptonff_dyturbo {
namespace {
std::mt19937_64 generator;
double uniform(){return ((generator()>>12)+.5)/4503599627370496.;}
}
std::string PhotosDecay::version(){return std::to_string(Photospp::Photos::VER_MAJOR)+"."+std::to_string(Photospp::Photos::VER_MINOR);}
PhotosDecay::PhotosDecay(double alpha,double cutoff,double precision):cutoff_(cutoff),precision_(precision) {
  if(!std::isfinite(alpha) || alpha<0 || !std::isfinite(cutoff) || cutoff<=0 || cutoff>=1 ||
     !std::isfinite(precision) || precision<=0 || precision>=1)throw std::runtime_error("Invalid PHOTOS numerical settings");
  using namespace Photospp;
  Photos::initialize();Photos::setRandomGenerator(uniform);
  Photos::setAlphaQED(alpha);Photos::setExponentiation(true);
  Photos::setDoubleBrem(false);Photos::setQuatroBrem(false);
  Photos::setInfraredCutOff(cutoff);Photos::setInterference(true);
  // Universal decay record: full production-dependent ME weights require
  // incoming flavours/spin information which this record does not supply.
  Photos::setMeCorrectionWtForW(false);Photos::setMeCorrectionWtForZ(false);
  Photos::setCorrectionWtForW(true);Photos::forceMass(13,.1056583755);
  phokey_.expeps=precision_;
}
std::pair<leptonff::FourVector,leptonff::FourVector> PhotosDecay::draw(
    leptonff::SpinDecay kind,const leptonff::BornDecayPoint& point,std::uint64_t seed) const {
  using namespace leptonff;using namespace Photospp;
  const bool neutral=kind==SpinDecay::neutral_vector;
  const auto& p=point.first;const auto& r=point.second;
  const FourVector q{p.energy+r.energy,p.px+r.px,p.py+r.py,p.pz+r.pz};
  const Real M=invariant_mass(p,r),m=.1056583755L;
  const auto rest=rest_frame_boost(p,q,-1);const Real norm=std::hypot(rest.pt(),rest.pz);
  const Real E=neutral?M/2:(M*M+m*m)/(2*M),k=std::sqrt(E*E-m*m);
  const FourVector a{E,k*rest.px/norm,k*rest.py/norm,k*rest.pz/norm};
  const FourVector b{M-E,-a.px,-a.py,-a.pz};
  PhotosHEPEVTEvent event;const int charge=kind==SpinDecay::wplus?1:-1;
  auto mother=new PhotosHEPEVTParticle(neutral?23:24*charge,2,0,0,0,M,M,-1,-1,-1,-1);
  auto first=new PhotosHEPEVTParticle(neutral?13:-13*charge,1,a.px,a.py,a.pz,a.energy,m,-1,-1,-1,-1);
  auto second=new PhotosHEPEVTParticle(neutral?-13:14*charge,1,b.px,b.py,b.pz,b.energy,neutral?m:0,-1,-1,-1,-1);
  event.addParticle(mother);event.addParticle(first);event.addParticle(second);
  first->setMothers({mother});second->setMothers({mother});mother->setDaughters({first,second});
  Photos::setInfraredCutOff(cutoff_);phokey_.expeps=precision_;
  generator.seed(seed);Photos::processParticle(mother);
  if(PH_HEPEVT_Interface::ME_channel!=0)throw std::runtime_error("Unexpected PHOTOS production-dependent ME channel");
  FourVector sum{};for(auto* item:event.getParticleList())if(item->getStatus()==1) {
    sum.energy+=item->getE();sum.px+=item->getPx();sum.py+=item->getPy();sum.pz+=item->getPz();
  }
  if(std::max({std::abs(sum.energy-M),std::abs(sum.px),std::abs(sum.py),std::abs(sum.pz)})>1e-7L*M)
    throw std::runtime_error("PHOTOS momentum conservation failure");
  const auto read=[&](PhotosHEPEVTParticle* v){return rest_frame_boost({v->getE(),v->getPx(),v->getPy(),v->getPz()},q,1);};
  return {read(first),read(second)};
}
}
