#include "leptonff/spin_decay.hpp"
#include "leptonff/coupled_lab.hpp"
#include "hard_primitive.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff {
namespace {
using Complex = std::complex<Real>;
using Matrix = std::array<Complex,16>;

Matrix identity() {
  Matrix result{};
  for (int i = 0; i < 4; ++i) result[5*i] = 1;
  return result;
}
Matrix scale(Matrix m, Real factor) {
  for (auto& x : m) x *= factor;
  return m;
}
Matrix add(const Matrix& a, const Matrix& b) {
  Matrix result{};
  for (int i = 0; i < 16; ++i) result[i] = a[i]+b[i];
  return result;
}
Matrix product(const Matrix& a, const Matrix& b) {
  Matrix result{};
  for (int i = 0; i < 4; ++i)
    for (int k = 0; k < 4; ++k) {
      const auto x = a[4*i+k];
      if (x == Complex{}) continue;
      for (int j = 0; j < 4; ++j) result[4*i+j] += x*b[4*k+j];
    }
  return result;
}
Matrix adjoint(const Matrix& a) {
  Matrix result{};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) result[4*i+j] = std::conj(a[4*j+i]);
  return result;
}
const std::array<Matrix,4>& gamma() {
  static const auto matrices = [] {
    std::array<Matrix,4> g{};
    g[0][0] = g[0][5] = 1;
    g[0][10] = g[0][15] = -1;
    const std::array<std::array<Complex,4>,3> pauli{{
        {{0,1,1,0}},{{0,Complex{0,-1},Complex{0,1},0}},{{1,0,0,-1}}}};
    for (int a = 0; a < 3; ++a)
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
          g[a+1][4*i+j+2] = pauli[a][2*i+j];
          g[a+1][4*(i+2)+j] = -pauli[a][2*i+j];
        }
    return g;
  }();
  return matrices;
}
Matrix slash(const FourVector& p) {
  return add(add(scale(gamma()[0],p.energy),scale(gamma()[1],-p.px)),
             add(scale(gamma()[2],-p.py),scale(gamma()[3],-p.pz)));
}
Real dot(const FourVector& a, const FourVector& b) {
  return a.energy*b.energy-a.px*b.px-a.py*b.py-a.pz*b.pz;
}


Matrix left_projector() {
  auto g5=product(product(product(gamma()[0],gamma()[1]),gamma()[2]),gamma()[3]);
  for(auto& x:g5)x*=Complex{0,1};
  return scale(add(identity(),scale(g5,-1)),.5L);
}
FourVector sum(const FourVector& a,const FourVector& b) {
  return {a.energy+b.energy,a.px+b.px,a.py+b.py,a.pz+b.pz};
}
std::array<Real,3> direction(const FourVector& p) {
  const Real norm=std::hypot(std::hypot(p.px,p.py),p.pz);
  if(norm==0)throw std::invalid_argument("Undefined rest-frame direction");
  return {p.px/norm,p.py/norm,p.pz/norm};
}
FourVector lightlike(Real e,const std::array<Real,3>& n) {
  return {e,e*n[0],e*n[1],e*n[2]};
}
void unit(const std::array<Real,3>& n) {
  Real square=0;
  for(auto x:n) {if(!std::isfinite(x))throw std::invalid_argument("Nonfinite spin direction");square+=x*x;}
  if(std::abs(square-1)>1e-9L)throw std::invalid_argument("Spin direction must be unit length");
}
Matrix radiative_current(SpinDecay kind,const SpinRadiativePoint& p,
                         const FourVector& e,int i) {
  const Real d1=2*dot(p.first,p.photon),d2=2*dot(p.second,p.photon);
  const auto first=add(scale(identity(),2*dot(p.first,e)),product(slash(e),slash(p.photon)));
  auto result=scale(product(first,gamma()[i]),1/d1);
  if(kind==SpinDecay::neutral_vector) {
    const auto second=add(scale(identity(),2*dot(p.second,e)),product(slash(p.photon),slash(e)));
    return add(result,scale(product(gamma()[i],second),-1/d2));
  }
  // gamma^i = slash(eta) requires eta^i=-1 for a spatial index.
  const FourVector eta{i==0?Real{1}:Real{0},i==1?Real{-1}:Real{0},
                       i==2?Real{-1}:Real{0},i==3?Real{-1}:Real{0}};
  const auto parent=sum(sum(p.first,p.second),p.photon);
  // Reduced unitary-gauge W graph, using P.eta=0 and k.e=0. It also
  // retains the longitudinal-propagator contribution when m_l != 0.
  const auto w=add(add(scale(gamma()[i],dot(parent,e)),
                         scale(slash(e),dot(p.photon,eta))),
                         scale(slash(p.photon),-dot(eta,e)));
  result=add(result,scale(w,-1/dot(parent,p.photon)));
  return product(result,left_projector());
}
SpinTensor squared_tensor(SpinDecay kind,const SpinRadiativePoint& p,Real m,
                         const FourVector& e) {
  const auto first=add(slash(p.first),scale(identity(),m));
  const auto second=add(slash(p.second),scale(identity(),kind==SpinDecay::neutral_vector?-m:0));
  std::array<Matrix,3> a,bar;
  for(int i=0;i<3;++i) {
    a[i]=radiative_current(kind,p,e,i+1);
    bar[i]=product(product(gamma()[0],adjoint(a[i])),gamma()[0]);
  }
  SpinTensor result{};
  for(int i=0;i<3;++i)for(int j=i;j<3;++j) {
    const auto square=product(product(product(first,a[i]),second),bar[j]);
    for(int k=0;k<4;++k)result[i][j]+=square[5*k];
    if(kind==SpinDecay::wplus)result[i][j]=std::conj(result[i][j]);
    result[j][i]=std::conj(result[i][j]);
  }
  return result;
}
SpinRadiativePoint kinematics(SpinDecay kind,Real M,Real m,Real z,Real c) {
  const Real r=m*m/(M*M),threshold=kind==SpinDecay::neutral_vector?4*r:r;
  if(!std::isfinite(M) || M<=0 || !std::isfinite(m) || m<0 ||
     !std::isfinite(z) || z<=threshold || z>=1 || !std::isfinite(c) ||
     std::abs(c)>1 || (m==0 && std::abs(c)==1))
    throw std::invalid_argument("Invalid spin radiative point");
  const Real root=std::sqrt(z),pair=M*root;
  const Real e1=kind==SpinDecay::neutral_vector?pair/2:(M*M*z+m*m)/(2*pair);
  const Real e2=pair-e1,k=std::sqrt(std::max(Real{0},(e1-m)*(e1+m)));
  const Real g=(1+z)/(2*root),v=-(1-z)/(1+z),px=k*std::sqrt((1-c)*(1+c));
  return {{g*(e1+v*k*c),px,0,g*(k*c+v*e1)},
          {g*(e2-v*k*c),-px,0,g*(-k*c+v*e2)},
          {M*(1-z)/2,0,0,M*(1-z)/2},{}};
}
} // namespace

Real projected_decay_hard_interval(SpinDecay kind,Real lo,Real hi,Real M,Real mu) {
  for(Real x:{lo,hi,M,mu})if(!std::isfinite(x))throw std::invalid_argument("Non-finite hard interval");
  if(lo<0 || hi<0 || M<=0 || mu<=0)throw std::invalid_argument("Invalid hard interval");
  if(lo>=1 || hi<=lo)return 0;
  if(hi==1)throw std::invalid_argument("Infrared unsafe hard interval");
  const bool charged=kind!=SpinDecay::neutral_vector;
  const Real log=2*std::log(M/mu);
  return (hi>1?0:detail::projected_hard_primitive(hi,charged,log))
    -detail::projected_hard_primitive(lo,charged,log);
}

Real spin_trace(const SpinTensor& t) {
  return (t[0][0]+t[1][1]+t[2][2]).real();
}
Real spin_contract(const SpinTensor& rho,const SpinTensor& t) {
  // Contract Hermitian pairs together. Summing nine complex terms separately
  // loses imaginary cancellation when real and projected soft terms are large.
  for(const auto* a:{&rho,&t}) {
    Real norm=1;
    for(auto row:*a)for(auto x:row)norm=std::max(norm,std::abs(x));
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
      if(!std::isfinite((*a)[i][j].real()) || !std::isfinite((*a)[i][j].imag()) ||
         std::abs((*a)[i][j]-std::conj((*a)[j][i]))>1e-9L*norm)
        throw std::runtime_error("Non-Hermitian spin tensor");
  }
  Real value=0;
  for(int i=0;i<3;++i) {
    value+=rho[i][i].real()*t[i][i].real();
    for(int j=i+1;j<3;++j)value+=2*(rho[i][j]*t[j][i]).real();
  }
  if(!std::isfinite(value))throw std::runtime_error("Invalid spin contraction");
  return value;
}
SpinTensor spin_born_tensor(SpinDecay kind,const std::array<Real,3>& n) {
  unit(n);SpinTensor result{};
  const Real sign=kind==SpinDecay::wminus?-1:kind==SpinDecay::wplus?1:0;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)result[i][j]=1.5L*((i==j?1:0)-n[i]*n[j]);
  // Sign follows the Dirac trace u(l-) gamma^i P_L v(nubar).
  result[0][1]+=Complex{0,sign*1.5L*n[2]};result[1][0]=std::conj(result[0][1]);
  result[1][2]+=Complex{0,sign*1.5L*n[0]};result[2][1]=std::conj(result[1][2]);
  result[2][0]+=Complex{0,sign*1.5L*n[1]};result[0][2]=std::conj(result[2][0]);
  return result;
}
SpinTensor spin_angular_projector(SpinDecay kind,const std::array<Real,3>& n) {
  const auto born=spin_born_tensor(kind,n);SpinTensor k{};
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    k[i][j]=Complex{(i==j?2:0)-5*n[i]*n[j],2*born[i][j].imag()/3};
  return k;
}
SpatialRotation decay_rotation(Real c,Real phi,Real twist) {
  if(!std::isfinite(c) || std::abs(c)>1 || !std::isfinite(phi) || !std::isfinite(twist))
    throw std::invalid_argument("Invalid decay rotation");
  const Real s=std::sqrt((1-c)*(1+c)),a=std::cos(phi),b=std::sin(phi),u=std::cos(twist),v=std::sin(twist);
  return {{{a*c*u-b*v,-a*c*v-b*u,a*s},
           {b*c*u+a*v,-b*c*v+a*u,b*s},{-s*u,s*v,c}}};
}
FourVector rotate_momentum(const FourVector& p,const SpatialRotation& r) {
  const std::array<Real,3> v{p.px,p.py,p.pz};std::array<Real,3> q{};
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)q[i]+=r[i][j]*v[j];
  return {p.energy,q[0],q[1],q[2]};
}
SpinTensor rotate_spin_tensor(const SpinTensor& t,const SpatialRotation& r) {
  SpinTensor q{};
  for(int i=0;i<3;++i)for(int j=i;j<3;++j) {
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)q[i][j]+=r[i][a]*t[a][b]*r[j][b];
    if(i==j)q[i][j]=q[i][j].real();
    q[j][i]=std::conj(q[i][j]);
  }
  return q;
}
FourVector rest_frame_boost(const FourVector& p,const FourVector& q,int sign) {
  const Real m2=dot(q,q);
  if(m2<=0 || q.energy<=0 || (sign!=1 && sign!=-1))throw std::invalid_argument("Invalid rest boost");
  const Real M=std::sqrt(m2),g=q.energy/M,bx=q.px/q.energy,by=q.py/q.energy,bz=q.pz/q.energy;
  const Real b2=bx*bx+by*by+bz*bz,bp=bx*p.px+by*p.py+bz*p.pz;
  const Real term=(b2>0?g*g*bp/(g+1):0)+sign*g*p.energy;
  return {g*(p.energy+sign*bp),p.px+term*bx,p.py+term*by,p.pz+term*bz};
}
SpinRadiativePoint spin_radiative_point(SpinDecay kind,Real M,Real m,Real z,Real c) {
  auto p=kinematics(kind,M,m,z,c);
  for(const auto e:{FourVector{0,1,0,0},FourVector{0,0,1,0}}) {
    const auto t=squared_tensor(kind,p,m,e);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)p.real[i][j]+=t[i][j];
  }
  const Real beta=kind==SpinDecay::neutral_vector?std::sqrt(1-4*m*m/(M*M*z)):1-m*m/(M*M*z);
  const Real factor=3*(1-z)*beta/(kind==SpinDecay::neutral_vector?16:8);
  for(auto& row:p.real)for(auto& x:row)x*=factor;
  return p;
}
Real spin_photon_ward_residual(SpinDecay kind,Real M,Real m,Real z,Real c) {
  const auto p=kinematics(kind,M,m,z,c);
  // Tests the actual spin-summed amplitude: the W graph cancels the charged
  // lepton graph before squaring. Neutral current cancels as a matrix.
  Real norm=0;
  for(int i=1;i<=3;++i)for(const auto x:radiative_current(kind,p,p.photon,i))norm=std::max(norm,std::abs(x));
  return norm;
}
Real projected_decay_hard_kernel(SpinDecay kind,Real z,Real M,Real mu) {
  if(!std::isfinite(z) || z<=0 || z>=1 || !std::isfinite(M) || M<=0 || !std::isfinite(mu) || mu<=0)
    throw std::invalid_argument("Invalid projected hard kernel");
  Real result=(1+z*z)/(1-z)*(2*std::log(M/mu)+std::log(z)+2*std::log1p(-z));
  if(kind!=SpinDecay::neutral_vector)result-=2*(1+z+z*z)/(3*(1-z));
  return result;
}
Real w_projected_decay_hard(const FourVector& mu,const FourVector& nu,Real scale,const WLabCuts& cuts) {
  cuts.validate();
  if(mu.pt()<=cuts.muon_pt_min || std::abs(mu.eta())>=cuts.muon_eta_max ||
     (cuts.neutrino_pt_min>0 && nu.pt()<=cuts.neutrino_pt_min))return 0;
  Real lo=cuts.muon_pt_min/mu.pt();
  if(cuts.transverse_mass_min>0)lo=std::max(lo,std::pow(cuts.transverse_mass_min/transverse_mass(mu,nu),2));
  return projected_decay_hard_interval(SpinDecay::wminus,lo,2,invariant_mass(mu,nu),scale);
}
Real z_projected_decay_hard(const FourVector& a,const FourVector& b,Real scale,const FiducialCuts& cuts) {
  cuts.validate();
  if(a.pt()<=cuts.pt_min || b.pt()<=cuts.pt_min || std::abs(a.eta())>=cuts.eta_max || std::abs(b.eta())>=cuts.eta_max)return 0;
  const Real M=invariant_mass(a,b),lo=std::pow(cuts.mass_min/M,2),hi=std::pow(cuts.mass_max/M,2);
  return projected_decay_hard_interval(SpinDecay::neutral_vector,std::max(cuts.pt_min/a.pt(),lo),hi,M,scale)
        +projected_decay_hard_interval(SpinDecay::neutral_vector,std::max(cuts.pt_min/b.pt(),lo),hi,M,scale);
}
std::array<SpinTensor,2> decay_hard_geometry_pair(SpinDecay kind,Real M,Real z,Real c,const SpatialRotation& rotation,
    const FourVector& parent,const DecayMeasurement& measurement,const DecayMeasurement& evolved,Real reference) {
  if(!std::isfinite(reference) || std::abs(dot(parent,parent)-M*M)>1e-8L*M*M)
    throw std::invalid_argument("Inconsistent decay parent or reference");
  const auto p=spin_radiative_point(kind,M,0,z,c);
  const auto real=rotate_spin_tensor(p.real,rotation);
  const auto lab=[&](const FourVector& v){return rest_frame_boost(rotate_momentum(v,rotation),parent,1);};
  const auto first_lab=lab(p.first),second_lab=lab(p.second);
  const Real f=measurement(first_lab,second_lab),fu=evolved(first_lab,second_lab);
  if(!std::isfinite(f) || !std::isfinite(fu))throw std::invalid_argument("Nonfinite decay measurement");
  std::array<SpinTensor,2> result{real,real};
  for(int k=0;k<2;++k)for(auto& row:result[k])for(auto& x:row)x*=k==0?f-reference:fu-f;
  const Real density=spin_trace(p.real)/3;
  const int sectors=kind==SpinDecay::neutral_vector?2:1;
  for(int sector=0;sector<sectors;++sector) {
    // Project the radiating cluster onto a massless two-body Born direction.
    auto n=direction(sector==0?sum(p.first,p.photon):p.first);
    const auto born=rotate_spin_tensor(spin_born_tensor(kind,n),rotation);
    auto first=lightlike(M/2,n);FourVector second{M/2,-first.px,-first.py,-first.pz};
    if(sector==0)first=first.scaled(z);else second=second.scaled(z);
    const Real weight=sectors==1?1:(sector==0?(1+c)/2:(1-c)/2);
    const auto a=lab(first),b=lab(second);
    const Real projected=measurement(a,b),projected_u=evolved(a,b);
    if(!std::isfinite(projected) || !std::isfinite(projected_u))throw std::invalid_argument("Nonfinite collinear measurement");
    for(int k=0;k<2;++k)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
      result[k][i][j]-=density*weight*born[i][j]*(k==0?projected-reference:projected_u-projected);
  }
  return result;
}
SpinTensor decay_hard_geometry(SpinDecay kind,Real M,Real z,Real c,const SpatialRotation& rotation,
    const FourVector& parent,const DecayMeasurement& measurement,Real reference) {
  return decay_hard_geometry_pair(kind,M,z,c,rotation,parent,measurement,measurement,reference)[0];
}
CoupledGeometry decay_hard_geometry_channels(SpinDecay kind,Real M,Real z,Real c,
    const SpatialRotation& rotation,const FourVector& parent,const ParentMeasurements& measurement,Real reference) {
  return decay_hard_geometry_channels_batch(kind,M,z,c,rotation,parent,{measurement},{reference})[0];
}
std::vector<CoupledGeometry> decay_hard_geometry_channels_batch(SpinDecay kind,Real M,Real z,Real c,
    const SpatialRotation& rotation,const FourVector& parent,
    const std::vector<ParentMeasurements>& measurements,const std::vector<Real>& references) {
  std::vector<std::array<Real,3>> controls;
  for(Real reference:references)controls.push_back({reference,reference,0});
  return decay_hard_geometry_linear_batch(kind,M,z,c,rotation,parent,measurements,controls);
}
std::vector<CoupledGeometry> decay_hard_geometry_linear_batch(SpinDecay kind,Real M,Real z,Real c,
    const SpatialRotation& rotation,const FourVector& parent,
    const std::vector<ParentMeasurements>& measurements,const std::vector<std::array<Real,3>>& references) {
  if(measurements.size()!=references.size())throw std::invalid_argument("Inconsistent geometry batch");
  for(const auto& reference:references)for(Real value:reference)
    if(!std::isfinite(value))throw std::invalid_argument("Nonfinite geometry reference");
  if(std::abs(dot(parent,parent)-M*M)>1e-8L*M*M)
    throw std::invalid_argument("Inconsistent coupled decay parent");
  const auto p=spin_radiative_point(kind,M,0,z,c);
  const auto real=rotate_spin_tensor(p.real,rotation);
  const auto lab=[&](const FourVector& v){return rest_frame_boost(rotate_momentum(v,rotation),parent,1);};
  const auto a=lab(p.first),b=lab(p.second),photon=lab(p.photon);
  std::vector<CoupledGeometry> out(measurements.size());
  for(std::size_t n=0;n<out.size();++n) {
    const auto f=measurements[n](a,b,photon);
    out[n]={real,real,real};
    for(int k=0;k<3;++k) {
      if(!std::isfinite(f[k]))throw std::invalid_argument("Nonfinite coupled measurement");
      for(auto& row:out[n][k])for(auto& x:row)x*=f[k]-references[n][k];
    }
  }
  const Real density=spin_trace(p.real)/3;
  const int sectors=kind==SpinDecay::neutral_vector?2:1;
  for(int sector=0;sector<sectors;++sector) {
    const auto n=direction(sector==0?sum(p.first,p.photon):p.first);
    const auto born=rotate_spin_tensor(spin_born_tensor(kind,n),rotation);
    auto first=lightlike(M/2,n);FourVector second{M/2,-first.px,-first.py,-first.pz};
    const auto photon=(sector==0?first:second).scaled(1-z);
    if(sector==0)first=first.scaled(z);else second=second.scaled(z);
    const auto fa=lab(first),fb=lab(second),fg=lab(photon);
    const Real weight=sectors==1?1:(sector==0?(1+c)/2:(1-c)/2);
    for(std::size_t n=0;n<out.size();++n) {
      const auto collinear=measurements[n](fa,fb,fg);
      for(int k=0;k<3;++k) {
        if(!std::isfinite(collinear[k]))throw std::invalid_argument("Nonfinite coupled collinear measurement");
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
          out[n][k][i][j]-=density*weight*born[i][j]*(collinear[k]-references[n][k]);
      }
    }
  }
  return out;
}
} // namespace leptonff
