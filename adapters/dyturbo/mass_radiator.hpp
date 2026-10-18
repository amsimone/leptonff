// SPDX-License-Identifier: MIT
#pragma once
#include "leptonff/photonic_dressing.hpp"
#include <algorithm>
#include <numeric>
#include <memory>
#include <stdexcept>
namespace leptonff_dyturbo {
enum class MassMatching { single_insertion, quadratic_finite };
// Inclusive mass reduction of the public two-leg operator. FFT products are
// zero padded: one global finite insertion, no circular large-loss wrap.
class MassRadiator {
  using Real=leptonff::Real;
  struct Table {
    Real h;std::vector<Real> w,prefix;
    Table(Real spacing,std::vector<Real> weights):h(spacing),w(std::move(weights)),prefix(w.size()) {
      for(std::size_t k=1;k<w.size();++k)prefix[k]=prefix[k-1]+w[k];
    }
    Real cumulative(Real v)const {
      if(v<=h/2)return 0;
      const auto k=std::size_t(v/h+.5L);
      if(k>=w.size())return prefix.back();
      return prefix[k-1]+w[k]*(v/h-k+.5L);
    }
    Real interval(Real lo,Real hi)const {
      if(lo>=1 || hi<=lo)return 0;
      if(lo<=0)throw std::invalid_argument("Mass radiator needs a positive observed lower bound");
      return cumulative(-std::log(lo))-cumulative(std::max(Real{0},-std::log(hi)))+(hi>1?w[0]:0);
    }
    // Moment of the SAME piecewise-constant log-loss density used by interval,
    // including its endpoint atom; not the atomic approximation exp(-s*k*h).
    Real moment(Real s)const {
      const Real x=s*h,cell=x==0?1:-std::expm1(-x)/x;
      Real value=w[0],correction=0;
      for(std::size_t k=1;k<w.size();++k) {
        const Real y=w[k]*std::exp(-s*(k-.5L)*h)*cell-correction,next=value+y;
        correction=(next-value)-y;value=next;
      }
      return value;
    }
  };
  static std::vector<Real> convolution(std::vector<Real> a,std::vector<Real> b) {
    const auto n=a.size();a.resize(2*n);b.resize(2*n);
    auto c=leptonff::circular_lattice_convolution(a,b).weights;c.resize(n);return c;
  }
  Real alpha_,scale_,L_;
  std::unique_ptr<Table> base_,scale_log_,scale_log2_;
public:
  MassRadiator(const leptonff::QEDParameters& q,const leptonff::LLGrid& grid,
               leptonff::PhotonicAccuracy accuracy=leptonff::PhotonicAccuracy::nll,
               MassMatching matching=MassMatching::single_insertion):
      alpha_(q.expansion_parameter()),scale_(q.mu_frag),L_(2*std::log(q.mu_frag/q.muon_mass)) {
    using namespace leptonff;
    const PhotonicLabDressing d(q,grid,1e-5L,accuracy);
    const auto& u=d.evolution().evolution();
    auto base=convolution(u.weights(),u.weights());
    const auto evolved_pair=matching==MassMatching::quadratic_finite?base:std::vector<Real>{};
    const auto insertion=convolution(d.nll_evolution()?d.nll_evolution()->finite_insertion_weights():
                                    d.evolution().evolved_boundary_weights(),u.weights());
    const auto hard=convolution(d.evolved_hard(),u.weights());
    const auto pp=convolution(convolution(u.generator_weights(),u.weights()),u.weights());
    auto logs=pp;
    for(std::size_t k=0;k<base.size();++k) {
      base[k]+=2*alpha_*(insertion[k]+hard[k]-d.evolution().first_order_weights()[k]-d.hard_source()[k]);
      logs[k]=2*alpha_*(pp[k]-u.generator_weights()[k]);
    }
    if(matching==MassMatching::quadratic_finite) {
      // Named NNLL completion: U^2 * a^2*C_pair^2/2 with
      // C_pair=2*(d1+h1+log(M^2/mu^2)*P0). No P1 cross terms.
      // The retained alpha^2 L^2 and alpha^2 L coefficients are unchanged.
      auto finite=d.evolution().initial_boundary_weights();
      for(std::size_t k=0;k<finite.size();++k)finite[k]+=d.hard_source()[k];
      const auto cc=convolution(convolution(finite,finite),evolved_pair);
      const auto pc=convolution(finite,pp);
      auto p2=convolution(u.generator_weights(),pp);
      for(std::size_t k=0;k<base.size();++k) {
        base[k]+=2*alpha_*alpha_*cc[k];
        logs[k]+=4*alpha_*alpha_*pc[k];
        p2[k]*=2*alpha_*alpha_;
      }
      scale_log2_=std::make_unique<Table>(grid.spacing(),std::move(p2));
    } else if(matching!=MassMatching::single_insertion)
      throw std::invalid_argument("Unknown mass matching prescription");
    base_=std::make_unique<Table>(grid.spacing(),std::move(base));
    scale_log_=std::make_unique<Table>(grid.spacing(),std::move(logs));
  }
  std::array<Real,3> bin(Real M,Real lower,Real upper)const {
    const auto lo=std::pow(lower/M,2),hi=std::pow(upper/M,2);
    if(lo>=1)return {};
    if(hi==1)throw std::invalid_argument("Strict Born endpoint veto");
    // Analytic whole-plus D1 integral, identical continuum convention to
    // LabFragmentation, independent of the averaged total-loss cells.
    const auto primitive=[&](Real x) {
      const Real t=1-x,l=std::log1p(-x);
      return (L_-1)*(-2*l-x-x*x/2)-2*(1.75L-(l*l-2*(t*l-t)+t*t*(l/2-.25L)));
    };
    const Real d1=(hi>1?0:primitive(hi))-primitive(lo);
    const Real hard=leptonff::projected_decay_hard_interval(leptonff::SpinDecay::neutral_vector,lo,hi,M,scale_);
    const Real correction=2*alpha_*(d1+hard),born=hi>1?1:0;
    const Real ell=2*std::log(M/scale_);
    return {born,born+correction,base_->interval(lo,hi)+ell*scale_log_->interval(lo,hi)+
            (scale_log2_?ell*ell*scale_log2_->interval(lo,hi):0)+correction};
  }
  std::array<Real,3> moment(Real M,Real N)const {
    if(!std::isfinite(M) || M<=0 || !std::isfinite(N) || N<1)
      throw std::invalid_argument("Mass radiator moment needs positive M and N>=1");
    const Real physical_log=L_+2*std::log(M/scale_);
    const leptonff::Distribution r1(0,{},{
      [physical_log](Real z) {return (1+z*z)/(1-z)*(physical_log+std::log(z)-1);}});
    const Real correction=2*alpha_*r1.moment(N);
    const Real ell=2*std::log(M/scale_);
    return {1,1+correction,base_->moment(N-1)+ell*scale_log_->moment(N-1)+
            (scale_log2_?ell*ell*scale_log2_->moment(N-1):0)+correction};
  }
};
}
