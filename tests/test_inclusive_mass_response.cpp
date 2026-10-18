#include "leptonff/photonic_dressing.hpp"
#include "check.hpp"
#include "../adapters/dyturbo/mass_radiator.hpp"
int main() {
  using namespace leptonff;
  // The mass-only angular remainder has zero trace at each radiation point.
  const FourVector parent{100,0,0,0};
  const auto bare=[](const FourVector& a,const FourVector& b) {
    const auto m=invariant_mass(a,b);return m>70 && m<95?Real{1}:Real{0};
  };
  const auto evolved=[](const FourVector& a,const FourVector& b) {
    return std::pow(invariant_mass(a,b)/100,2);
  };
  for(Real z:{.2L,.6L,.9L})for(Real c:{-.7L,.3L}) {
    const auto tensors=decay_hard_geometry_pair(SpinDecay::neutral_vector,100,z,c,
      decay_rotation(.3L,1.2L,.7L),parent,bare,evolved);
    for(const auto& tensor:tensors)close(spin_trace(tensor),0,1e-10L,
      "inclusive mass radiator has no angular-remainder trace");
  }
  LLGrid grid;grid.bins=32768;grid.log_span=32;grid.soft_log_cutoff=grid.spacing()/2;
  QEDParameters q;q.alpha=1/137.035999084L;q.mu_frag=91.1876L;
  const leptonff_dyturbo::MassRadiator fast(q,grid);
  const PhotonicLabDressing reference(q,grid,1e-5L,PhotonicAccuracy::nll);
  const leptonff_dyturbo::MassRadiator fast_ll(q,grid,PhotonicAccuracy::ll);
  const PhotonicLabDressing reference_ll(q,grid,1e-5L,PhotonicAccuracy::ll);
  const Real M=91.123L,a=q.expansion_parameter();
  close(fast.moment(M,1)[2],1,2e-11L,"complete radiator number moment");
  close(fast.moment(M,2)[1],1+a*(44.L/9-8.L/3*2*std::log(M/q.muon_mass)),2e-10L,
    "complete radiator first mass moment matches independent decay integral");
  throws<std::invalid_argument>([&]{fast.moment(M,.5L);},"reject unsupported mass moment");
  const leptonff_dyturbo::MassRadiator quadratic(q,grid,PhotonicAccuracy::nll,
    leptonff_dyturbo::MassMatching::quadratic_finite);
  close(quadratic.moment(M,1)[2],1,3e-11L,"quadratic completion preserves normalization");
  close(quadratic.moment(M,2)[1],fast.moment(M,2)[1],0,"quadratic completion preserves first order");
  const FourVector p{M/2,M/2,0,0},r{M/2,-M/2,0,0};
  for(const auto interval:std::array<std::pair<Real,Real>,3>{{{60,75},{80,90},{90,95}}}) {
    const auto [lo,hi]=interval;const auto reduced=fast.bin(M,lo,hi);
    const FiducialCuts cuts{M*1e-4L/2,100,lo,hi};
    const auto ff=reference.fragmentation().z(p,r,cuts);
    const Real hard=z_projected_decay_hard(p,r,q.mu_frag,cuts);
    close(reduced[0],ff.born,1e-12L,"mass reduction Born");
    close(reduced[1],ff.fixed_order+a*hard,1e-10L,"mass reduction continuum first order");
    close(reduced[2],ff.matched+a*(hard+reference.z_projected_shift(p,r,cuts)),2e-4L,
      "total-loss mass reduction converges to public two-leg cell operator");
    const auto ll=reference_ll.fragmentation().z(p,r,cuts);
    close(fast_ll.bin(M,lo,hi)[2],ll.matched+a*(hard+reference_ll.z_projected_shift(p,r,cuts)),2e-4L,
      "LL mass reduction preserves the global finite-insertion convention");
  }
  // Vary evolution and hard matching together at a fixed physical parent mass.
  for(Real factor:{.5L,2.L}) {
    auto varied=q;varied.mu_frag*=factor;
    for(auto accuracy:{PhotonicAccuracy::ll,PhotonicAccuracy::nll}) {
      const leptonff_dyturbo::MassRadiator shifted(varied,grid,accuracy);
      const PhotonicLabDressing public_shifted(varied,grid,1e-5L,accuracy);
      for(const auto interval:std::array<std::pair<Real,Real>,3>{{{60,75},{80,90},{90,95}}}) {
        const auto [lo,hi]=interval;const auto v=shifted.bin(M,lo,hi);
        close(v[1],fast.bin(M,lo,hi)[1],2e-14L,"fragmentation scale cancels at order alpha");
        const FiducialCuts cuts{M*1e-4L/2,100,lo,hi};
        const auto ff=public_shifted.fragmentation().z(p,r,cuts);
        const Real h=z_projected_decay_hard(p,r,varied.mu_frag,cuts);
        close(v[2],ff.matched+a*(h+public_shifted.z_projected_shift(p,r,cuts)),2e-4L,
          "varied-scale mass reduction agrees with public matching");
      }
    }
  }
  q.alpha=0;const leptonff_dyturbo::MassRadiator zero(q,grid);
  close(zero.moment(M,5)[2],1,1e-12L,"zero-alpha complete moment");
  for(const auto interval:std::array<std::pair<Real,Real>,2>{{{80,90},{90,95}}}) {
    const auto v=zero.bin(M,interval.first,interval.second);
    close(v[1],v[0],1e-12L,"mass reduction zero-alpha first order");
    close(v[2],v[0],1e-12L,"mass reduction zero-alpha resummation");
  }
  q.alpha=1e-4L;const leptonff_dyturbo::MassRadiator small(q,grid);
  q.alpha/=2;const leptonff_dyturbo::MassRadiator half(q,grid);
  const auto s=small.bin(M,80,90),t=half.bin(M,80,90);
  close((s[2]-s[1])/(t[2]-t[1]),4,.04L,
    "inclusive mass reduction has no unmatched linear term");
}
