#include "leptonff/vector_decay.hpp"
#include "leptonff/vector_current.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const Real pi = std::acos(-1.0L);
  close(real_dilogarithm(0),0,1e-15L,"Li2 zero");
  close(real_dilogarithm(1),pi*pi/6,1e-14L,"Li2 endpoint");
  close(real_dilogarithm(-1),-pi*pi/12,1e-14L,"Li2 negative endpoint");
  close(real_dilogarithm(.5L),pi*pi/12-std::pow(std::log(2.0L),2)/2,
        1e-14L,"Li2 half");
  for (Real x : {-100.0L,-.25L,.2L,.9L})
    close(real_dilogarithm(x),integrate([=](Real u) {
      return -std::log1p(-x*u)/u;
    },0,1).value,1e-10L,"Li2 vs independent integral");

  QEDParameters q;
  const Real M = 91.1876L;
  const auto hard = vector_mass_hard_coefficient(M,M);
  close(hard.moment(1),1.5L,1e-12L,"inclusive finite hard coefficient");
  close(hard.moment(2),61.0L/6,1e-8L,"independent hard N2 moment");
  const Real L = 2*std::log(M/q.muon_mass);
  close(vector_mass_leading_power_coefficient(M,q.muon_mass).moment(2),
        115.0L/18-8*L/3,1e-8L,"independent full spectrum N2 moment");

  for (const MassWindow window : {MassWindow{60,120},MassWindow{70,85},MassWindow{0,120}}) {
    Real central = 0;
    for (Real scale : {M/2,M,2*M}) {
      q.mu_frag = scale;
      const auto result = vector_mass_window(M,window,q);
      close(result.matched_coefficient,result.leading_power_coefficient,2e-8L,
            "finite hard plus FF reproduces physical leading-power spectrum");
      if (scale == M) central = result.leading_power_rate;
      if (scale == 2*M)
        close(result.leading_power_rate,central,1e-10L,"fragmentation scale cancellation");
    }
    // The exact massive expression is independently normalised to massless Born.
    // Allow logarithmically enhanced power corrections rather than demanding Cm^2.
    for (Real muon : {.1056583755L,.25L,.5L,1.0L,2.0L}) {
      q.muon_mass = muon;
      q.mu_0 = muon;
      q.mu_frag = M;
      const auto result = vector_mass_window(M,window,q);
      const Real r = std::pow(muon/M,2);
      close(result.massive_rate,result.leading_power_rate,
            100*q.expansion_parameter()*r*(1+std::abs(std::log(r)))+10*r*r,
            "massive minus leading-power rate is power suppressed");
    }
  }
  q.muon_mass = .1056583755L;
  q.mu_frag = M;
  Real starting_scale_reference = 0;
  for (Real factor : {.5L,1.0L,2.0L}) {
    q.mu_0 = factor*q.muon_mass;
    const auto value = vector_mass_window(M,{60,120},q).matched_coefficient;
    if (factor == .5L) starting_scale_reference = value;
    close(value,starting_scale_reference,1e-9L,"hard+FF starting-scale cancellation");
  }
  close(massive_vector_inclusive_coefficient(M,1e-6L),1.5L,1e-10L,
        "massive inclusive high-energy limit");
  close(massive_vector_born_ratio(M,1e-6L),1,1e-14L,"massive Born high-energy limit");
  close(massive_vector_real_cumulative(M,.1L,0),0,1e-15L,"below pair threshold");
  for (Real muon : {.1056583755L,1.0L,2.0L,25.0L})
    for (Real z : {.45L,.7L,.93L}) {
      const auto point = vector_radiative_point(M,muon,z,.3L);
      close(point.minus.energy+point.plus.energy+point.photon.energy,M,1e-12L,
            "radiative energy conservation");
      close(point.minus.pz+point.plus.pz+point.photon.pz,0,1e-12L,
            "radiative longitudinal conservation");
      close(point.minus.energy*point.minus.energy-point.minus.px*point.minus.px
            -point.minus.pz*point.minus.pz,muon*muon,2e-12L,"massive muon shell");
      close(invariant_mass(point.minus,point.plus),M*std::sqrt(z),1e-12L,
            "radiative dimuon mass");
      close(vector_current_ward_residual(M,muon,z,.3L),0,1e-12L,"photon Ward identity");
      close(vector_real_angular_integral(M,muon,z,{1e-7L,1e-8L,8192}),
            massive_vector_real_spectrum(M,muon,z),2e-6L,
            "Dirac trace vs independent exact massive spectrum");
    }
  if (!(massive_vector_real_cumulative(M,.1L,.8L)>
        massive_vector_real_cumulative(M,.1L,.5L)))
    throw std::runtime_error("Massive real cumulative not monotonic");
  throws<std::invalid_argument>([&] { MassWindow{60,M}.validate(M); },
                                "infrared unsafe Born endpoint cut");
  throws<std::invalid_argument>([&] { massive_vector_real_cumulative(M,.1L,1); },
                                "divergent real cumulative endpoint");
  throws<std::invalid_argument>([&] { massive_vector_born_ratio(M,M/2); },
                                "pair threshold");
  throws<std::invalid_argument>([] { real_dilogarithm(1.1L); },"complex Li2 branch");
  std::cout << "Vector mass-spectrum hard/FF matching and massive checks passed\n";
}
