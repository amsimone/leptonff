#include "leptonff/drell_yan.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const SpatialRotation identity{{{1,0,0},{0,1,0},{0,0,1}}};
  DrellYanAngularWeights w;w.unpolarised=-2;w.angular={.4L,.1L,-.2L,.3L,-.1L,.07L,-.03L,.09L};
  for(auto kind:{SpinDecay::wplus,SpinDecay::wminus,SpinDecay::neutral_vector}) {
    const auto rho=production_spin_tensor(kind,w,identity);
    close(spin_trace(rho),w.unpolarised,1e-15L,"production tensor trace is unnormalised rate");
    for(Real c:{-.8L,.1L,.7L})for(Real phi:{.2L,1.3L,2.4L}) {
      const Real s=std::sqrt(1-c*c);const std::array<Real,3> n{s*std::cos(phi),s*std::sin(phi),c};const auto& h=w.angular;
      Real polynomial=w.unpolarised*(1+c*c)+h[0]*(1-3*c*c)/2+h[1]*2*s*c*std::cos(phi)+h[2]*s*s*std::cos(2*phi)/2
        +h[5]*s*s*std::sin(2*phi)+h[6]*2*s*c*std::sin(phi);
      if(kind!=SpinDecay::neutral_vector)polynomial+=h[3]*s*std::cos(phi)+h[4]*c+h[7]*s*std::sin(phi);
      close(spin_contract(rho,spin_born_tensor(kind,n)),.75L*polynomial,2e-14L,"independent angular polynomial including W odd channels");
    }
    w.unpolarised=0;const auto zero=production_spin_tensor(kind,w,identity);
    close(spin_trace(zero),0,1e-15L,"zero normalisation does not divide by Born");w.unpolarised=-2;
  }
  for(const FourVector q:std::vector<FourVector>{{100,0,0,0},{125,30,20,60},{150,0,0,100}}) {
    const auto r=collins_soper_rotation(q);const auto rho=production_spin_tensor(SpinDecay::wplus,w,r);
    close(spin_trace(rho),w.unpolarised,2e-14L,"recoil rotation preserves production trace");
  }
  auto bad=identity;bad[0][0]=2;
  throws<std::invalid_argument>([&]{production_spin_tensor(SpinDecay::wplus,w,bad);},"reject malformed angular frame");
  std::cout<<"Unnormalised Drell-Yan angular weights, signed/zero terms and explicit recoil frame passed\n";
}
