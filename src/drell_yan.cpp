#include "leptonff/drell_yan.hpp"
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
using Axis=std::array<Real,3>;
Axis unit(Axis x) {const Real n=std::hypot(std::hypot(x[0],x[1]),x[2]);if(n==0 || !std::isfinite(n))throw std::invalid_argument("Degenerate angular axis");for(auto& v:x)v/=n;return x;}
Real dot(const Axis& a,const Axis& b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
}
SpinTensor production_spin_tensor(SpinDecay kind,const DrellYanAngularWeights& weights,const SpatialRotation& rotation) {
  if(!std::isfinite(weights.unpolarised))throw std::invalid_argument("Non-finite production coefficient");
  for(Real x:weights.angular)if(!std::isfinite(x))throw std::invalid_argument("Non-finite angular coefficient");
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    Real norm=0;for(int k=0;k<3;++k) {if(!std::isfinite(rotation[k][i]))throw std::invalid_argument("Non-finite angular frame");norm+=rotation[k][i]*rotation[k][j];}
    if(std::abs(norm-(i==j?1:0))>1e-10L)throw std::invalid_argument("Angular frame is not orthonormal");
  }
  const Real determinant=rotation[0][0]*(rotation[1][1]*rotation[2][2]-rotation[1][2]*rotation[2][1])
    -rotation[0][1]*(rotation[1][0]*rotation[2][2]-rotation[1][2]*rotation[2][0])
    +rotation[0][2]*(rotation[1][0]*rotation[2][1]-rotation[1][1]*rotation[2][0]);
  if(std::abs(determinant-1)>1e-10L)throw std::invalid_argument("Angular frame must be right handed");
  const auto& h=weights.angular;const Real U=weights.unpolarised;
  SpinTensor rho{};rho[0][0]=(2*U-h[0]-h[2])/4;rho[1][1]=(2*U-h[0]+h[2])/4;rho[2][2]=h[0]/2;
  rho[0][1]=-h[5]/2;rho[0][2]=-h[1]/2;rho[1][2]=-h[6]/2;
  if(kind!=SpinDecay::neutral_vector) {
    const Real sign=spin_born_tensor(kind,{0,0,1})[0][1].imag()/1.5L;
    rho[0][1]+=std::complex<Real>{0,h[4]/(4*sign)};
    rho[0][2]+=std::complex<Real>{0,-h[7]/(4*sign)};
    rho[1][2]+=std::complex<Real>{0,h[3]/(4*sign)};
  }
  for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)rho[j][i]=std::conj(rho[i][j]);
  return rotate_spin_tensor(rho,rotation);
}
SpatialRotation collins_soper_rotation(const FourVector& q) {
  if(!std::isfinite(q.energy) || q.energy<=0 || q.energy*q.energy-q.px*q.px-q.py*q.py-q.pz*q.pz<=0)
    throw std::invalid_argument("Timelike boson required for angular frame");
  const auto p1=rest_frame_boost({1,0,0,1},q,-1),p2=rest_frame_boost({1,0,0,-1},q,-1);
  const auto a=unit({p1.px,p1.py,p1.pz}),b=unit({p2.px,p2.py,p2.pz});
  const auto z=unit({a[0]-b[0],a[1]-b[1],a[2]-b[2]});Axis x{a[0]+b[0],a[1]+b[1],a[2]+b[2]};
  if(dot(x,x)<1e-24L) {const Axis reference=std::abs(z[0])<.9L?Axis{1,0,0}:Axis{0,1,0};for(int i=0;i<3;++i)x[i]=reference[i]-dot(reference,z)*z[i];}
  x=unit(x);const auto y=unit({z[1]*x[2]-z[2]*x[1],z[2]*x[0]-z[0]*x[2],z[0]*x[1]-z[1]*x[0]});
  SpatialRotation r{};for(int i=0;i<3;++i){r[i][0]=x[i];r[i][1]=y[i];r[i][2]=z[i];}return r;
}
} // namespace leptonff
