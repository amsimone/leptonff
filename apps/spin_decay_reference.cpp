#include "leptonff/spin_decay.hpp"
#include "leptonff/fragmentation.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace leptonff;
namespace {
constexpr Real M=80.385L,pt=25,eta=2.4L;
Real acceptance(Real momentum) {
  return momentum>pt?std::min(std::tanh(eta),std::sqrt(1-pt*pt/(momentum*momentum))):0;
}
// Independent rotational reduction of the unpolarised W measurement. The
// massive calculation subtracts the Born measurement and needs NO collinear
// cutoff; the massless geometry subtracts the z-projected Born measurement.
Real radiation(Real mass,Real soft,Real angular) {
  const Real r=mass*mass/(M*M),born=acceptance(M*(1-r)/2);
  return integrate([&](Real z) {
    const Real projected=mass>0?born:acceptance(z*M/2),beta=mass>0?(z-r)/(z+r):1;
    const Real lower=mass>0?-std::atanh(beta):-1+angular;
    const Real upper=mass>0?std::atanh(beta):1-angular;
    std::vector<Real> breaks;
    for(Real momentum:{pt,pt*std::cosh(eta)}) {
      const Real energy=std::sqrt(momentum*momentum+mass*mass);
      const Real c=((1+z)*(z+r)-4*z*energy/M)/((1-z)*(z-r));
      if(c>-1 && c<1)breaks.push_back(mass>0?std::atanh(beta*c):c);
    }
    return integrate([&](Real u) {
      const Real c=mass>0?std::tanh(u)/beta:u;
      const auto p=spin_radiative_point(SpinDecay::wminus,M,mass,z,c);
      const Real momentum=std::sqrt(std::max(Real{0},p.first.energy*p.first.energy-mass*mass));
      const Real jacobian=mass>0?1/(beta*std::pow(std::cosh(u),2)):1;
      return spin_trace(p.real)/3*(acceptance(momentum)-projected)*jacobian;
    },lower,upper,{1e-5L,1e-5L,8192},breaks).value;
  },std::max(r,1e-7L),1-soft,{1e-4L,1e-4L,8192},{2*pt/M}).value;
}
}
int main() {
  try {
    QEDParameters q;q.mu_frag=M;
    const Real born=acceptance(M/2);
    const auto measurement=[](Real z){return acceptance(z*M/2);};
    const Real ff=fragmentation_coefficient_order1(Species::muon_minus,Species::muon_minus,q)
      .apply(measurement,0,{}, {2*pt/M});
    const Real projected=integrate([&](Real z) {
      return projected_decay_hard_kernel(SpinDecay::wminus,z,M,q.mu_frag)*(measurement(z)-born);
    },0,1,{1e-7L,1e-7L,8192},{2*pt/M}).value;
    const Real coarse=radiation(0,1e-5L,1e-5L),fine=radiation(0,1e-6L,1e-6L);
    const Real massive=radiation(q.muon_mass,1e-6L,0),full=ff+projected+fine;
    const Real r=std::pow(q.muon_mass/M,2),massive_born=std::pow(1-r,2)*(1+r/2)*acceptance(M*(1-r)/2);
    if(std::abs(fine-coarse)>5e-4L || std::abs(full-massive)>.003L)
      throw std::runtime_error("W massive/leading-power matching or cutoff refinement failed");
    std::cout<<std::setprecision(17)<<"{\n  \"version\": \""<<LEPTONFF_VERSION<<"\",\n"
      <<"  \"scope\": \"unitary unpolarised W rest decay; exact massive real versus LP hard+FF\",\n"
      <<"  \"parent_mass\": "<<M<<",\n  \"muon_mass\": "<<q.muon_mass<<",\n"
      <<"  \"pt_min\": "<<pt<<",\n  \"eta_max\": "<<eta<<",\n"
      <<"  \"born\": "<<born<<",\n  \"ff_coefficient1\": "<<ff<<",\n"
      <<"  \"projected_hard_coefficient1\": "<<projected<<",\n"
      <<"  \"geometry_coarse\": "<<coarse<<",\n  \"geometry_refined\": "<<fine<<",\n"
      <<"  \"geometry_refinement_shift\": "<<fine-coarse<<",\n"
      <<"  \"full_coefficient1\": "<<full<<",\n  \"massive_coefficient1\": "<<massive<<",\n"
      <<"  \"LP_minus_massive_coefficient\": "<<full-massive<<",\n"
      <<"  \"LP_rate\": "<<born+q.expansion_parameter()*full<<",\n"
      <<"  \"massive_rate\": "<<massive_born+q.expansion_parameter()*massive<<",\n"
      <<"  \"integration_abs_rel_inner\": 1e-5,\n  \"integration_abs_rel_outer\": 1e-4\n}\n";
  } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
