#include "leptonff/moment_evolution.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  QEDFlavourScheme one;one.fermions={{"mu",-1,.1056583755L,1,0}};
  const std::vector<bool> active{true};
  for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)for(Real N:{2.L,3.L,5.L}) {
    close(flavour_leading_kernel(one,active,i,j).moment(N),leading_kernel(one_family_species[i],one_family_species[j]).moment(N),2e-11L,"one-family LO kernel reduction");
    close(flavour_next_to_leading_kernel(one,active,i,j).moment(N),next_to_leading_kernel(one_family_species[i],one_family_species[j]).moment(N),2e-9L,"one-family NLO kernel reduction");
  }
  auto charged=one;charged.fermions.push_back({"up",2.L/3,1.5L,3,0});charged.fermions.push_back({"down",-1.L/3,4.75L,3,0});
  const std::vector<bool> all(3,true);const auto sums=qed_charge_sums(charged,all);
  close(sums.squared,8.L/3,1e-15L,"closed Dirac loop charge sum includes colours");
  close(sums.fourth,44.L/27,1e-15L,"two-loop charge fourth powers");
  for(const auto& scheme:{one,charged,QEDFlavourScheme::leptons()}) {
    const std::vector<bool> mask(scheme.fermions.size(),true);
    for(std::size_t f=0;f<scheme.fermions.size();++f)
      close(flavour_valence_kernel(scheme,mask,f).moment(1),0,1e-9L,"all-flavour valence charge conservation at N1");
    for(std::size_t i=0;i<scheme.size();++i) {
      Real lo=0,nlo=0,boundary=0;
      for(std::size_t j=0;j<scheme.size();++j) {
        lo+=flavour_leading_kernel(scheme,mask,i,j).moment(2);
        nlo+=flavour_next_to_leading_kernel(scheme,mask,i,j).moment(2);
        boundary+=flavour_initial_condition(scheme,mask,i,j,10).moment(2);
      }
      close(lo,0,1e-10L,"all-flavour per-parent LO momentum sum");
      close(nlo,0,1e-8L,"all-flavour per-parent NLO momentum sum");
      close(boundary,0,1e-8L,"all-flavour boundary momentum sum");
    }
  }
  for(Real z:{.07L,.3L,.8L}) {
    const Real singlet=nlo_lepton_component(NLOLeptonComponent::singlet_pair).resolved_density(z);
    const Real forward=flavour_next_to_leading_kernel(charged,all,0,2).resolved_density(z);
    const Real reverse=flavour_next_to_leading_kernel(charged,all,2,0).resolved_density(z);
    close(forward,4*singlet/3,1e-12L,"produced-colour multiplicity, not closed-loop nf");
    close(forward,3*reverse,1e-12L,"parent versus daughter colour factor");
    close(flavour_next_to_leading_kernel(charged,all,0,3).resolved_density(z),forward,1e-12L,"unlike-flavour pair charge symmetry");
  }
  // Source Eq. (27), hep-ph/0504192: the finite time-like heavy component
  // does NOT vanish when matching at the mass. This also checks its log.
  for(Real factor:{.5L,1.L,2.L})for(Real N:{2.L,3.L,5.L}) {
    const auto& heavy=charged.fermions[1];const Real mu=factor*heavy.mass,L=2*std::log(factor),q2=4.L/9;
    const Real p=(N*N+N+2)/(N*(N*N-1));
    const Real analytic=q2*(p*(L-1)+4/std::pow(N-1,2)-4/(N*N)+2/std::pow(N+1,2));
    close(flavour_threshold_coefficient(charged,1,2,charged.photon(),mu).moment(N),analytic,2e-10L,"time-like threshold analytic Mellin coefficient");
    close(flavour_threshold_coefficient(charged,1,charged.photon(),charged.photon(),mu).moment(N),-2*3*q2*L/3,1e-14L,"photon decoupling colour/charge and sign");
  }
  if(flavour_threshold_coefficient(one,0,0,one.photon(),one.fermions[0].mass).moment(2)==0)
    throw std::runtime_error("PDF-like zero-at-mass threshold negative control");
  const auto leptons=QEDFlavourScheme::leptons();const Real alpha=1/137.035999084L,mu=.1056583755L;
  const Real converted=leptonic_alpha_msbar(alpha,mu,leptons);
  close(1/converted,1/alpha-2*std::log(mu/leptons.fermions[0].mass)/(3*std::acos(-1.L)),1e-12L,"one-loop OS/MSbar input matching");
  throws<std::invalid_argument>([&]{leptonic_alpha_msbar(alpha,mu,charged);},"no fabricated on-shell quark/hadronic conversion");
  throws<std::invalid_argument>([]{auto s=QEDFlavourScheme::leptons();s.fermions[0].charge=0;s.validate();},"invalid charged flavour");
  throws<std::invalid_argument>([]{auto s=QEDFlavourScheme::leptons();s.fermions[2].matching_scale=100;s.validate();},"matching scale near mass enforced");
  std::cout<<"Flavour charge/colour factors, momentum sums, time-like matching and input conversion passed\n";
}
