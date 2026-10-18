// Optional PHOTOS decay comparison on the SAME signed MCFM contributions.
// PHOTOS has global state: this adapter deliberately requires one host thread.
#include "Photos.h"
#include "PhotosHEPEVTEvent.h"
#include "PhotosHEPEVTParticle.h"
#include "photos_comparison.hpp"
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
namespace {
using namespace leptonff;
using namespace Photospp;
std::mt19937_64 rng;
double random_number() {return ((rng()>>11)+.5)/9007199254740992.;}
Real maximum_conservation_error=0;
double linear_scale=1;
bool first_order=false;
struct Init {
  Init() {
    const char* threads=std::getenv("OMP_NUM_THREADS");
    if(!threads || std::string(threads)!="1")
      throw std::invalid_argument("PHOTOS comparison requires OMP_NUM_THREADS=1");
    Photos::initialize();Photos::setRandomGenerator(random_number);
    const char* alpha=std::getenv("LEPTONFF_ALPHA");
    const char* order=std::getenv("LEPTONFF_PHOTOS_ORDER");
    const bool exponentiation=!(order && std::string(order)=="1");
    first_order=!exponentiation;
    const char* proposal=std::getenv("LEPTONFF_PHOTOS_LINEAR_SCALE");
    linear_scale=proposal?std::stod(proposal):1;
    if(!std::isfinite(linear_scale) || linear_scale<1 || (exponentiation && linear_scale!=1))
      throw std::invalid_argument("PHOTOS linear scale must be >=1 and requires first order");
    Photos::setAlphaQED((alpha?std::stod(alpha):1/137.035999084)/linear_scale);
    Photos::setExponentiation(exponentiation);Photos::setDoubleBrem(false);Photos::setQuatroBrem(false);
    const char* cut=std::getenv("LEPTONFF_PHOTOS_CUTOFF");
    const double cutoff=cut?std::stod(cut):1e-5;
    Photos::setInfraredCutOff(cutoff);Photos::setInterference(true);
    Photos::setMeCorrectionWtForW(true);Photos::setMeCorrectionWtForZ(false);
    Photos::forceMass(13,.1056583755);
    std::cerr<<"leptonff PHOTOS 3.56 comparison: exponentiation "<<(exponentiation?"on":"off (first order)")
      <<", W ME on, Z ME off, interference on; common LEPTONFF_ALPHA, cutoff="<<cutoff
      <<", linear proposal scale="<<linear_scale<<"; one thread\n";
  }
  ~Init(){std::cerr<<"leptonff PHOTOS max rest-frame momentum error GeV="<<maximum_conservation_error<<'\n';}
};
FourVector boost(const FourVector& p,const FourVector& q,int sign) {
  const Real mass=std::sqrt(q.energy*q.energy-q.px*q.px-q.py*q.py-q.pz*q.pz);
  const Real gamma=q.energy/mass;
  const Real bx=q.px/q.energy,by=q.py/q.energy,bz=q.pz/q.energy;
  const Real b2=bx*bx+by*by+bz*bz,dot=bx*p.px+by*p.py+bz*p.pz;
  const Real term=(b2>0?(gamma-1)*dot/b2:0)+sign*gamma*p.energy;
  return {gamma*(p.energy+sign*dot),p.px+term*bx,p.py+term*by,p.pz+term*bz};
}
std::uint64_t seed(const double* a,const double* b) {
  std::uint64_t result=1469598103934665603ULL;
  for(const auto* p:{a,b})for(int k=0;k<4;++k) {
    std::uint64_t bits;std::memcpy(&bits,p+k,sizeof(bits));
    result^=bits;result*=1099511628211ULL;
  }
  return result;
}
bool selected(int process,const FourVector& a,const FourVector& b) {
  if(a.pt()<=25 || std::abs(a.eta())>=2.4)return false;
  if(process!=31)return true;
  const Real mass=invariant_mass(a,b);
  return b.pt()>25 && std::abs(b.eta())<2.4 && mass>60 && mass<120;
}
}

std::vector<std::array<leptonff::Real,4>> leptonff_photos_batch(int process,const double* p3,const double* p4,
    const double* ref3,const double* ref4,const std::vector<leptonff::DecayMeasurement>& measurements) {
  static Init init;
  using namespace leptonff;
  using namespace Photospp;
  const auto from=[](const double* p){return FourVector{p[3],p[0],p[1],p[2]};};
  auto first=from(process==1?p4:p3),second=from(process==1?p3:p4);
  // The hard host supplies massless leptons; preserve three-momenta when
  // repairing its tiny mapping residuals, as in the analytic FF adapter.
  first.energy=std::hypot(first.pt(),first.pz);second.energy=std::hypot(second.pt(),second.pz);
  const FourVector q{first.energy+second.energy,first.px+second.px,
      first.py+second.py,first.pz+second.pz};
  const Real mass=invariant_mass(first,second),m=.1056583755L;
  if(mass<=2*m)throw std::runtime_error("PHOTOS parent below muon threshold");
  const auto rest=boost(first,q,-1);
  const Real norm=std::hypot(rest.pt(),rest.pz);
  const Real energy=process==31?mass/2:(mass*mass+m*m)/(2*mass);
  const Real k=std::sqrt(energy*energy-m*m);
  const FourVector a{energy,k*rest.px/norm,k*rest.py/norm,k*rest.pz/norm};
  const FourVector b{mass-energy,-a.px,-a.py,-a.pz};
  std::vector<std::array<Real,4>> result(measurements.size());
  const auto ba=boost(a,q,1),bb=boost(b,q,1);
  for(std::size_t i=0;i<result.size();++i)result[i][0]=measurements[i](ba,bb);
  PhotosHEPEVTEvent event;
  const int charge=process==1?1:-1;
  auto mother=new PhotosHEPEVTParticle(process==31?23:24*charge,2,0,0,0,mass,mass,-1,-1,-1,-1);
  auto lepton=new PhotosHEPEVTParticle(process==31?13:-charge*13,1,a.px,a.py,a.pz,a.energy,m,-1,-1,-1,-1);
  auto other=new PhotosHEPEVTParticle(process==31?-13:charge*14,1,b.px,b.py,b.pz,b.energy,process==31?m:0,-1,-1,-1,-1);
  event.addParticle(mother);event.addParticle(lepton);event.addParticle(other);
  lepton->setMothers({mother});other->setMothers({mother});mother->setDaughters({lepton,other});
  // Real and ALL mapped dipoles share the reference real momenta and therefore
  // exactly the same random stream. No independent photon draw on counterevents.
  rng.seed(seed(ref3,ref4));
  Photos::processParticle(mother);
  const auto read=[](PhotosHEPEVTParticle* p){return FourVector{p->getE(),p->getPx(),p->getPy(),p->getPz()};};
  const auto ra=boost(read(lepton),q,1),rb=boost(read(other),q,1);
  // With one photon, W ME on / Z ME off, the transition is linear in alpha.
  // Sampling at alpha/s avoids PHOTOS's crude-probability bound at high mass;
  // F_B+s*(F_after-F_B) is an unbiased estimator at the requested alpha.
  for(std::size_t i=0;i<result.size();++i) {
    const Real raw=measurements[i](ra,rb);
    result[i][1]=result[i][0]+linear_scale*(raw-result[i][0]);
    result[i][2]=result[i][1]-result[i][0];
  }
  FourVector sum{0,0,0,0};int photons=0;
  for(auto* p:event.getParticleList())if(p->getStatus()==1) {
    sum.energy+=p->getE();sum.px+=p->getPx();sum.py+=p->getPy();sum.pz+=p->getPz();
    if(p->getPdgID()==22)++photons;
  }
  const Real error=std::max({std::abs(sum.energy-mass),std::abs(sum.px),std::abs(sum.py),std::abs(sum.pz)});
  maximum_conservation_error=std::max(maximum_conservation_error,error);
  if(first_order && photons>1)throw std::runtime_error("PHOTOS first-order mode emitted multiple photons");
  if(error>1e-7L*mass)throw std::runtime_error("PHOTOS momentum conservation failed");
  for(auto& r:result)r[3]=photons>0?1:0;
  return result;
}

void leptonff_photos_comparison(int process,const double* p3,const double* p4,
    const double* ref3,const double* ref4,double* result) {
  const auto values=leptonff_photos_batch(process,p3,p4,ref3,ref4,
    {[=](const leptonff::FourVector& a,const leptonff::FourVector& b)->leptonff::Real {
      return selected(process,a,b)?1:0;
    }});
  for(int k=0;k<4;++k)result[k]=double(values[0][k]);
}
