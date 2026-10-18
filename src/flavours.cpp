#include "leptonff/flavours.hpp"
#include "flavour_kernels.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace leptonff {
void QEDFlavourScheme::validate() const {
  if(fermions.empty() || fermions.size()>16 || muon_flavour>=fermions.size())
    throw std::invalid_argument("Invalid charged-flavour basis");
  std::set<std::string> names;
  for(const auto& f:fermions) {
    if(f.name.empty() || !names.insert(f.name).second || !std::isfinite(f.charge) ||
       f.charge==0 || std::abs(f.charge)>2 || f.colours<1 || f.colours>3 ||
       !std::isfinite(f.mass) || !std::isfinite(f.matching_scale) ||
       (f.pointlike?(f.colours!=3 || f.mass!=0 || f.matching_scale<1):
        (f.mass<=0 || f.matching_scale<0 || f.threshold()/f.mass<.25L || f.threshold()/f.mass>4)))
      throw std::invalid_argument("Invalid fermion mass/charge/matching scale");
  }
  const auto& mu=fermions[muon_flavour];
  if(mu.charge!=-1 || mu.colours!=1)throw std::invalid_argument("Observed muon must have charge -1 and one colour");
}
std::vector<bool> QEDFlavourScheme::active(Real scale) const {
  validate();if(!std::isfinite(scale) || scale<=0)throw std::invalid_argument("Invalid flavour scale");
  std::vector<bool> result;for(const auto& f:fermions)result.push_back(scale>=f.threshold());return result;
}
QEDFlavourScheme QEDFlavourScheme::leptons(Real factor) {
  if(!std::isfinite(factor) || factor<.25L || factor>4)throw std::invalid_argument("Invalid heavy threshold factor");
  QEDFlavourScheme s;
  s.fermions={{"e",-1,.00051099895069L,1,0},{"mu",-1,.1056583755L,1,0},
              {"tau",-1,1.77693L,1,1.77693L*factor}};
  s.muon_flavour=1;s.validate();return s;
}
QEDFlavourScheme QEDFlavourScheme::hadronic_pointlike(Real scale) {
  if(!std::isfinite(scale) || scale<1 || scale>4)throw std::invalid_argument("Hadronic matching scale must be 1..4 GeV");
  auto s=leptons();s.include_gluon=true;
  for(auto f:std::vector<ChargedFermion>{{"u",2.L/3,0,3,scale,true},{"d",-1.L/3,0,3,scale,true},
       {"s",-1.L/3,0,3,scale,true},{"c",2.L/3,1.67L,3,std::max(scale,1.67L)},
       {"b",-1.L/3,4.78L,3,4.78L}})s.fermions.push_back(f);
  s.validate();return s;
}
QEDFlavourScheme QEDFlavourScheme::from_file(const std::string& path) {
  std::ifstream input(path);if(!input)throw std::invalid_argument("Cannot read flavour file: "+path);
  QEDFlavourScheme s;std::string line;bool muon=false;
  while(std::getline(input,line)) {
    line=line.substr(0,line.find('#'));std::istringstream row(line);ChargedFermion f;
    if(!(row>>f.name))continue;
    if(!(row>>f.charge>>f.colours>>f.mass>>f.matching_scale))throw std::invalid_argument("Incomplete flavour-file row");
    std::string extra;if(row>>extra){if(extra!="pointlike" || (row>>extra))throw std::invalid_argument("Extra flavour-file field");f.pointlike=true;}
    if(f.name=="mu"){s.muon_flavour=s.fermions.size();muon=true;}
    s.fermions.push_back(f);
  }
  if(!muon)throw std::invalid_argument("Flavour file must identify the observed muon as mu");
  s.validate();return s;
}
QEDChargeSums qed_charge_sums(const QEDFlavourScheme& s,const std::vector<bool>& active) {
  s.validate();if(active.size()!=s.fermions.size())throw std::invalid_argument("Invalid active-flavour mask");
  QEDChargeSums sums;
  for(std::size_t i=0;i<active.size();++i)if(active[i]) {
    const auto& f=s.fermions[i];const Real q2=f.charge*f.charge;
    sums.squared+=f.colours*q2;sums.fourth+=f.colours*q2*q2;
  }
  return sums;
}
namespace detail {
Distribution flavour_basis_distribution(std::size_t k) {
  using S=Species;using C=NLOLeptonComponent;using G=NLOPhotonComponent;
  const QEDParameters boundary{0,1,1,1};
  switch(k) {
    case Pff:return leading_kernel(S::muon_minus,S::muon_minus);
    case Pfg:return leading_kernel(S::muon_minus,S::photon);
    case Pgf:return leading_kernel(S::photon,S::muon_minus);
    case Delta:return {1};
    case Photonic:return nlo_lepton_component(C::photonic);
    case SameInterference:return nlo_lepton_component(C::same_charge_interference);
    case OppositeInterference:return nlo_lepton_component(C::opposite_charge_interference);
    case NonSingletPair:return nlo_lepton_component(C::nonsinglet_pair);
    case SingletPair:return nlo_lepton_component(C::singlet_pair);
    case Pfg1:return nlo_photon_component(G::fermion_to_photon);
    case Pgf1Direct:return nlo_photon_component(G::photon_to_fermion_direct);
    case Pgf1Vacuum:return nlo_photon_component(G::photon_to_fermion_vacuum);
    case Pgg1:return nlo_photon_component(G::photon_to_photon);
    case QCDgg:return {-0.5L,[](Real z){return 6*((1-z)/z+z*(1-z));},
      {[](Real z){return 6*z/(1-z);}},1};
    case Dff0:return initial_condition_coefficient(S::muon_minus,S::muon_minus,boundary);
    case Dfg0:return initial_condition_coefficient(S::muon_minus,S::photon,boundary);
  }
  throw std::invalid_argument("Invalid flavour kernel basis element");
}
namespace {
void indices(const QEDFlavourScheme& s,std::size_t i,std::size_t j) {
  if(i>=s.size() || j>=s.size())throw std::invalid_argument("Flavour channel index out of range");
}
}
FlavourExpansion flavour_kernel_expansion(const QEDFlavourScheme& s,const std::vector<bool>& active,
    std::size_t i,std::size_t j,int order) {
  indices(s,i,j);const auto sums=qed_charge_sums(s,active);FlavourExpansion c{};
  if(s.include_gluon && (i==s.gluon() || j==s.gluon()))return c;
  const auto g=s.photon();
  if((i!=g && !active[i/2]) || (j!=g && !active[j/2]))return c;
  if(order!=0 && order!=1)throw std::invalid_argument("Invalid kernel order");
  const Real qi=i==g?0:s.fermions[i/2].charge*s.fermions[i/2].charge;
  const Real qj=j==g?0:s.fermions[j/2].charge*s.fermions[j/2].charge;
  if(order==0) {
    if(i==g && j==g)c[Delta]=-sums.beta0();
    else if(i==g)c[Pgf]=s.fermions[j/2].colours*qj;
    else if(j==g)c[Pfg]=qi;
    else if(i==j)c[Pff]=qi;
    return c;
  }
  if(i==g && j==g)c[Pgg1]=sums.fourth;
  else if(i==g) {
    c[Pgf1Direct]=s.fermions[j/2].colours*qj*qj;
    c[Pgf1Vacuum]=s.fermions[j/2].colours*qj*sums.squared;
  }else if(j==g)c[Pfg1]=qi*qi;
  else {
    // Produced-flavour colour factor belongs to the daughter; closed loops
    // multiply NS only. Identical-flavour interference has no flavour sum.
    c[SingletPair]=s.fermions[j/2].colours*qi*qj;
    if(i==j){c[Photonic]=qi*qi;c[SameInterference]=qi*qi;c[NonSingletPair]=qi*sums.squared;}
    else if(i/2==j/2)c[OppositeInterference]=qi*qi;
  }
  return c;
}
FlavourExpansion flavour_boundary_expansion(const QEDFlavourScheme& s,const std::vector<bool>& active,
    std::size_t i,std::size_t j,Real mu) {
  indices(s,i,j);qed_charge_sums(s,active);
  if(!std::isfinite(mu) || mu<=0)throw std::invalid_argument("Invalid boundary scale");
  const auto g=s.photon();FlavourExpansion c{};
  if(s.include_gluon && (i==s.gluon() || j==s.gluon()))return c;
  if((i!=g && !active[i/2]) || (j!=g && !active[j/2]))return c;
  if(i==g && j==g) {
    for(std::size_t f=0;f<active.size();++f)if(active[f])
      c[Delta]-=4*s.fermions[f].colours*std::pow(s.fermions[f].charge,2)*std::log(mu/(s.fermions[f].pointlike?s.fermions[f].threshold():s.fermions[f].mass))/3;
  }else if(i==g) {
    const auto& f=s.fermions[j/2];c[Pgf]=2*f.colours*f.charge*f.charge*std::log(mu/(f.pointlike?f.threshold():f.mass));
  }else if(i==j || j==g) {
    const auto& f=s.fermions[i/2];const Real q=f.charge*f.charge;
    if(!f.pointlike)c[i==j?Dff0:Dfg0]=q;
    c[i==j?Pff:Pfg]=2*q*std::log(mu/(f.pointlike?f.threshold():f.mass));
  }
  return c;
}
FlavourExpansion flavour_matching_expansion(const QEDFlavourScheme& s,std::size_t heavy,
    std::size_t i,std::size_t j,Real mu) {
  indices(s,i,j);s.validate();
  if(heavy>=s.fermions.size() || !std::isfinite(mu) || mu<=0)throw std::invalid_argument("Invalid threshold matching request");
  if(s.fermions[heavy].pointlike || (s.include_gluon && (i==s.gluon() || j==s.gluon())))return {};
  const auto& f=s.fermions[heavy];const Real q=f.charge*f.charge,L=2*std::log(mu/f.mass);
  FlavourExpansion c{};
  if(i==s.photon() && j==s.photon())c[Delta]=-2.L*f.colours*q*L/3;
  else if(i!=s.photon() && i/2==heavy && j==s.photon()){c[Dfg0]=q;c[Pfg]=q*L;}
  return c;
}
FlavourExpansion qcd_kernel_expansion(const QEDFlavourScheme& s,const std::vector<bool>& active,std::size_t i,std::size_t j) {
  indices(s,i,j);qed_charge_sums(s,active);FlavourExpansion c{};
  if(!s.include_gluon)return c;
  const auto g=s.gluon();const auto quark=[&](std::size_t k){return k<s.photon() && active[k/2] && s.fermions[k/2].colours==3;};
  if(i==g && j==g) {
    c[QCDgg]=1;for(std::size_t f=0;f<active.size();++f)if(active[f] && s.fermions[f].colours==3)c[Delta]-=1.L/3;
  }else if(i==g && quark(j))c[Pgf]=.5L;
  else if(quark(i) && j==g)c[Pfg]=4.L/3;
  else if(i==j && quark(i))c[Pff]=4.L/3;
  return c;
}
} // namespace detail
namespace {
Distribution distribution(const detail::FlavourExpansion& c) {
  Distribution d;for(std::size_t k=0;k<c.size();++k)if(c[k]!=0)d=d+detail::flavour_basis_distribution(k).scaled(c[k]);return d;
}
}
Distribution flavour_leading_kernel(const QEDFlavourScheme& s,const std::vector<bool>& a,std::size_t i,std::size_t j) {
  return distribution(detail::flavour_kernel_expansion(s,a,i,j,0));
}
Distribution flavour_next_to_leading_kernel(const QEDFlavourScheme& s,const std::vector<bool>& a,std::size_t i,std::size_t j) {
  return distribution(detail::flavour_kernel_expansion(s,a,i,j,1));
}
Distribution flavour_valence_kernel(const QEDFlavourScheme& s,const std::vector<bool>& a,std::size_t f) {
  const auto sums=qed_charge_sums(s,a);if(f>=s.fermions.size())throw std::invalid_argument("Invalid valence flavour");
  if(!a[f])return {};
  const Real q2=s.fermions[f].charge*s.fermions[f].charge;
  return (nlo_lepton_component(NLOLeptonComponent::photonic)
    +nlo_lepton_component(NLOLeptonComponent::same_charge_interference)
    +nlo_lepton_component(NLOLeptonComponent::opposite_charge_interference).scaled(-1)).scaled(q2*q2)
    +nlo_lepton_component(NLOLeptonComponent::nonsinglet_pair).scaled(q2*sums.squared);
}
Distribution flavour_initial_condition(const QEDFlavourScheme& s,const std::vector<bool>& a,std::size_t i,std::size_t j,Real mu) {
  return distribution(detail::flavour_boundary_expansion(s,a,i,j,mu));
}
Distribution flavour_threshold_coefficient(const QEDFlavourScheme& s,std::size_t h,std::size_t i,std::size_t j,Real mu) {
  return distribution(detail::flavour_matching_expansion(s,h,i,j,mu));
}
Distribution qcd_leading_kernel(const QEDFlavourScheme& s,const std::vector<bool>& a,std::size_t i,std::size_t j) {
  return distribution(detail::qcd_kernel_expansion(s,a,i,j));
}
Real leptonic_alpha_msbar(Real alpha,Real mu,const QEDFlavourScheme& s) {
  s.validate();if(!std::isfinite(alpha) || alpha<0 || alpha>=1)throw std::invalid_argument("Invalid on-shell alpha");
  Real c=0;const auto active=s.active(mu);
  for(std::size_t f=0;f<s.fermions.size();++f) {
    if(s.fermions[f].colours!=1)throw std::invalid_argument("On-shell conversion requires leptonic content; no hadronic model supplied");
    if(active[f])c+=4*std::pow(s.fermions[f].charge,2)*std::log(mu/s.fermions[f].mass)/3;
  }
  const Real denominator=1-alpha*c/(2*std::acos(-1.L));
  if(denominator<=0 || alpha/denominator>=1)throw std::invalid_argument("OS/MSbar conversion outside perturbative range");
  return alpha/denominator;
}
} // namespace leptonff
