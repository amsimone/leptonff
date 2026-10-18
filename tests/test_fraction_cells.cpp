#include "leptonff/lab_fragmentation.hpp"
#include "check.hpp"
#include <algorithm>
#include <numeric>

int main() {
  using namespace leptonff;
  const Real h=.3L,minimum=.001L;
  std::vector<Real> w(32,0),v(32,0);w[0]=.7L;w[1]=-.2L;w[2]=.4L;w[4]=-.1L;
  v[0]=-.3L;v[1]=.6L;v[3]=-.2L;v[5]=.1L;
  const SignedFractionTable t(h,w,minimum,FractionMeasure::conservative_cells);
  const auto literal_interval=[&](const SignedFractionTable& table,Real lo,Real hi,bool inverse) {
    Real out=lo<1 && hi>1?table.weights()[0]:0;
    for(std::size_t j=1;j<table.weights().size();++j) {
      const Real a=std::max(lo,table.cell_lower(j)),b=std::min(hi,table.cell_upper(j));
      if(b>a)out+=table.weights()[j]/(table.cell_upper(j)-table.cell_lower(j))*(inverse?std::log(b/a):b-a);
    }
    return out;
  };
  for(Real lo:{minimum,.15L,.5L,.8L,.99L,1.L})for(Real hi:{.3L,.6L,1.L,1.1L,2.L}) {
    close(t.interval(lo,hi),literal_interval(t,lo,hi,false),2e-13L,"signed cell interval vs literal cell overlap");
    close(t.inverse_interval(lo,hi),literal_interval(t,lo,hi,true),2e-13L,"inverse prefix vs literal logarithmic integral");
  }
  Real energy=0;std::vector<Real> edges;
  for(std::size_t j=0;j<t.fractions().size();++j) {
    energy+=t.weights()[j]*t.fractions()[j];
    if(j) {
      close((t.cell_lower(j)+t.cell_upper(j))/2,t.fractions()[j],1e-14L,"cell energy mean at original node");
      if(j>1)close(t.cell_upper(j),t.cell_lower(j-1),1e-14L,"resolved cells meet without gaps");
      if(t.cell_lower(j)>minimum)edges.push_back(t.cell_lower(j));
      if(j==1)edges.push_back(t.cell_upper(j));
    }
  }
  close(t.interval(minimum),std::accumulate(w.begin(),w.end(),Real{0}),1e-14L,"number conserved without renormalisation");
  const Real numeric_energy=minimum*t.interval(minimum)+integrate([&](Real z){return t.interval(z);},minimum,1,{},edges).value;
  close(numeric_energy,energy,2e-12L,"prefix measurement preserves signed energy moment");
  close(t.interval(.999L,1),0,0,"strict upper endpoint excludes exact atom");
  close(t.interval(.999L),w[0],1e-14L,"unresolved endpoint atom remains exact");
  const Real middle=t.fractions()[1],epsilon=1e-8L;
  close(t.interval(middle-epsilon)-t.interval(middle+epsilon),
        2*epsilon*w[1]/(t.cell_upper(1)-t.cell_lower(1)),1e-13L,"resolved CDF continuous across node");

  for(auto ma:{FractionMeasure::atoms,FractionMeasure::conservative_cells})
  for(auto mb:{FractionMeasure::atoms,FractionMeasure::conservative_cells}) {
    const SignedFractionTable a(h,w,minimum,ma),b(h,v,minimum,mb);
    const auto interval=[&](const SignedFractionTable& table,std::size_t i,Real lo,Real hi) {
      if(table.measure()==FractionMeasure::atoms || i==0)return lo<table.fractions()[i] && table.fractions()[i]<hi?1.L:0.L;
      return std::max(0.L,std::min(hi,table.cell_upper(i))-std::max(lo,table.cell_lower(i)))/(table.cell_upper(i)-table.cell_lower(i));
    };
    for(Real l1:{minimum,.4L,.8L})for(Real l2:{minimum,.3L,.9L})
    for(const auto window:{std::pair<Real,Real>{0,2},{.1L,.4L},{.35L,1},{.5L,1.1L}}) {
      const Real lo=window.first,hi=window.second;Real literal=0;
      for(std::size_t i=0;i<a.weights().size();++i)if(a.weights()[i])
      for(std::size_t j=0;j<b.weights().size();++j)if(b.weights()[j]) {
        const auto& wa=a.weights()[i];const auto& wb=b.weights()[j];
        if(ma==FractionMeasure::atoms || i==0) {
          const Real z=a.fractions()[i];if(z>l1)literal+=wa*wb*interval(b,j,std::max(l2,lo/z),hi/z);
        }else if(mb==FractionMeasure::atoms || j==0) {
          const Real z=b.fractions()[j];if(z>l2)literal+=wa*wb*interval(a,i,std::max(l1,lo/z),hi/z);
        }else {
          const Real al=std::max(l1,a.cell_lower(i)),au=a.cell_upper(i);
          const Real bl=std::max(l2,b.cell_lower(j)),bu=b.cell_upper(j);
          if(al>=au || bl>=bu)continue;
          std::vector<Real> contacts;
          for(Real q:{lo/bl,lo/bu,hi/bl,hi/bu})if(q>al+1e-12L && q<au-1e-12L)contacts.push_back(q);
          std::sort(contacts.begin(),contacts.end());
          contacts.erase(std::unique(contacts.begin(),contacts.end(),[](Real x,Real y){return std::abs(x-y)<1e-12L;}),contacts.end());
          const Real area=integrate([&](Real z){return std::max(0.L,std::min(bu,hi/z)-std::max(bl,lo/z));},al,au,{},contacts).value;
          literal+=wa*wb*area/((a.cell_upper(i)-a.cell_lower(i))*(b.cell_upper(j)-b.cell_lower(j)));
        }
      }
      close(a.pair(b,l1,l2,lo,hi),literal,3e-11L,"exact cell pair vs independent rectangle integration");
      close(a.pair(b,l1,l2,lo,hi),b.pair(a,l2,l1,lo,hi),3e-11L,"signed pair symmetry, both atom/cell orientations");
    }
  }
  // Independent continuum benchmark: D(z)=1. Project z*D with the analytic
  // linear-z basis of the evolution grid, then measure its sharp cuts.
  // Conservative cells must approach 1-cut at second order on this smooth
  // regular source, rather than merely reproducing a discrete cell formula.
  Real olderror=0;
  for(Real spacing:{.2L,.1L,.05L}) {
    const auto n=static_cast<std::size_t>(std::ceil(16/spacing))+1;
    std::vector<Real> weights(n,0);
    for(std::size_t j=0;j<n;++j) {
      const Real high=std::exp(-spacing*j),low=std::exp(-spacing*(j+1));
      weights[j]+=(high-low)*(2*high+low)/(6*high);
      if(j+1<n)weights[j+1]+=(high-low)*(high+2*low)/(6*low);
    }
    const SignedFractionTable cells(spacing,weights,1e-5L,FractionMeasure::conservative_cells);
    const SignedFractionTable atoms(spacing,weights,1e-5L);
    Real error=0,atomic_error=0;
    for(Real cut:{.21L,.37L,.53L,.68L,.79L}) {
      error=std::max(error,std::abs(cells.interval(cut)-(1-cut)));
      atomic_error=std::max(atomic_error,std::abs(atoms.interval(cut)-(1-cut)));
    }
    if(olderror)close(olderror/error,4,.08L,"second-order CDF convergence to independent continuum");
    if(error>=atomic_error)throw std::runtime_error("Cell continuum fixture did not improve the atomic cut error");
    olderror=error;
  }
  throws<std::invalid_argument>([&]{t.interval(minimum/2);},"cell support remains enforced");
  std::cout<<"Signed conservative fraction cells, exact endpoint, moments and hyperbolic pair cuts passed\n";
}
