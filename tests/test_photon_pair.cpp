#include "leptonff/photon_pair.hpp"
#include "leptonff/vector_current.hpp"
#include "leptonff/vector_decay.hpp"
#include "check.hpp"
#include <algorithm>

int main() {
  using namespace leptonff;
  const Real M = 91.1876L;
  const IntegrationOptions tight{1e-9L,1e-9L,8192};
  QEDParameters q;
  const auto H = vector_photon_pair_mass_hard_coefficient(M,M);
  close(H.moment(2,tight),-13.0L/3,1e-8L,"analytic photon-pair hard energy moment");
  throws<std::invalid_argument>([&] { H.moment(1); },"divergent photon number moment");
  // Integrate the independently traced massive current at fixed photon-muon
  // invariant. Its coordinates, density and numerical integral differ from
  // the analytic t primitive used by the new massive spectrum.
  for (Real m : {.1056583755L,2.0L}) for (Real y : {.2L,.6L,.85L}) {
    const Real r = std::pow(m/M,2), v = 1-y;
    const Real root = std::sqrt(v*v-4*r);
    const Real tmin = 2*r*y/(v-2*r+root);
    const Real tmax = y*(v-2*r+root)/(2*(y+r));
    const Real trace = integrate([&](Real logt) {
      const Real t = std::exp(logt), z = 1-y-t;
      const Real beta = std::sqrt(1-4*r/z);
      const Real c = std::clamp((y-t)/((y+t)*beta),Real{-1},Real{1});
      return vector_real_density(M,m,z,c)*2*t/((y+t)*beta);
    },std::log(tmin),std::log(tmax),{1e-7L,1e-8L,8192}).value;
    close(massive_vector_photon_pair_spectrum(M,m,y),trace,2e-6L,
          "photon-pair spectrum vs integrated Dirac current");
  }
  for (Real y : {.2L,.6L,.85L}) {
    Real previous = 1;
    for (Real factor : {1.0L,.5L,.25L}) {
      q.muon_mass = .1056583755L*factor;
      q.mu_0 = q.muon_mass;
      const Real lp = H.resolved_density(y)+fragmentation_coefficient_order1(
          Species::muon_minus,Species::photon,q).resolved_density(y);
      const Real error = std::abs(massive_vector_photon_pair_spectrum(M,q.muon_mass,y)-lp);
      if (error > .35L*previous || error > .005L)
        throw std::runtime_error("Photon-pair massive limit fails quadratic convergence");
      previous = error;
    }
  }
  for (Real m : {.1056583755L,2.0L}) {
    const Real pair_energy = integrate([&](Real y) {
      return y*massive_vector_photon_pair_spectrum(M,m,y);
    },0,1-2*m/M,{1e-7L,1e-8L,8192}).value;
    const Real photon_energy = integrate([&](Real z) {
      return (1-z)*massive_vector_real_spectrum(M,m,z);
    },4*std::pow(m/M,2),1,{1e-7L,1e-8L,8192}).value;
    close(2*pair_energy,photon_energy,2e-6L,
          "both specified-charge pair moments reproduce photon energy");
  }
  HardDecayOptions options;
  options.integration = {1e-7L,1e-7L,8192};
  options.collinear_cutoff = 1e-8L;
  for (auto fractions : {std::pair<Real,Real>{1,1},{.9L,.95L}}) {
    const Real fg = fractions.first, fm = fractions.second;
    const FiducialCuts cuts;
    for (Real y : {.55L,.7L,.9L})
      close(vector_photon_pair_acceptance(M,y,0,fg,fm,cuts),
          isotropic_acceptance(M,y*fg,fm,cuts),1e-14L,"photon local collinear mapping");
    const auto central = fiducial_vector_photon_pair_hard(M,M,cuts,fg,fm,options);
    auto refined = options;
    refined.collinear_cutoff /= 4;
    close(fiducial_vector_photon_pair_hard(M,M,cuts,fg,fm,refined).total,
        central.total,2e-6L,"finite photon geometry sliver refinement");
    q = QEDParameters{};
    const auto B = [&](Real y) { return isotropic_acceptance(M,y*fg,fm,cuts); };
    const Real low = std::pow(cuts.mass_min/M,2)/(fg*fm);
    const auto ff = [&](Real mu) {
      q.mu_frag = mu;
      return fragmentation_coefficient_order1(Species::muon_minus,Species::photon,q)
          .apply(B,low,tight,{2*cuts.pt_min/(M*fg),fm/fg});
    };
    const Real full = central.total+ff(M);
    for (Real scale : {.5L,2.0L})
      close(fiducial_vector_photon_pair_hard(M,scale*M,cuts,fg,fm,options).total
          +ff(scale*M),full,1e-7L,"photon-parent hard plus FF scale cancellation");
  }
  const FiducialCuts inclusive{0,100,60,120};
  const auto r = fiducial_vector_photon_pair_hard(M,M,inclusive,1,1,options);
  close(r.geometry,0,0,"angular-inclusive photon hard geometry vanishes");
  close(r.total,H.apply([](Real) { return Real{1}; },std::pow(60/M,2),tight),
      2e-7L,"photon mass-only reduction");
  close(fiducial_vector_photon_pair_hard(M,M,FiducialCuts{},0,1).total,0,0,
      "zero photon fraction");
  throws<std::invalid_argument>([&] {
    fiducial_vector_photon_pair_hard(M,M,{0,100,0,120});
  },"unresolved same-direction pair excluded");
  throws<std::invalid_argument>([&] {
    vector_photon_pair_acceptance(M,.8L,.3L,1,1,FiducialCuts{});
  },"unphysical photon-pair invariants");
  std::cout << "Photon-parent hard coefficient, massive current and scale checks passed\n";
}
