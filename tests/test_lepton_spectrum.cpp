#include "leptonff/lepton_spectrum.hpp"
#include "leptonff/hard_dressing.hpp"
#include "check.hpp"
#include <algorithm>

int main() {
  using namespace leptonff;
  const Real M=91.1876L;
  const FiducialCuts cuts;
  // One accepted muon can exceed the extra threshold while its partner does
  // not. Tightening BOTH cuts would incorrectly give zero here.
  const SingleLeptonSelection pt35{SingleLeptonCut::pt_above,35};
  close(isotropic_acceptance(M,90/M,60/M,cuts,pt35),
      std::sqrt(1-std::pow(25.0L/30,2))/2,1e-14L,"one-lepton mean, not pair-minimum pT");
  // Signed non-endpoint weights catch indexing errors hidden by Born-only
  // and telescoping histogram checks. Compare a literal unsymmetrised tensor.
  const std::vector<Real> a{1,-.2L,.3L},b{.8L,.4L,-.1L};
  const Real spacing=.15L;
  for(SingleLeptonSelection selection:{pt35,{SingleLeptonCut::absolute_eta_below,.5L}}) {
    Real direct=0;
    for(std::size_t i=0;i<a.size();++i) for(std::size_t j=0;j<b.size();++j)
      direct+=a[i]*b[j]*isotropic_acceptance(M,std::exp(-spacing*i),std::exp(-spacing*j),cuts,selection);
    close(decay_lattice_acceptance(M,cuts,spacing,a,b,false,selection),direct,1e-14L,
          "selected signed lattice vs literal non-endpoint tensor");
  }

  for(SingleLeptonSelection selection:{pt35,{SingleLeptonCut::absolute_eta_below,.5L}})
    for(Real f1:{.75L,1.0L}) {
      const Real f2=.9L;
      std::vector<Real> breaks;
      for(Real f:{f1,f2}) {
        const Real p=M*f/2;
        const Real fid=std::min(std::tanh(cuts.eta_max),std::sqrt(1-std::pow(cuts.pt_min/p,2)));
        breaks.push_back(fid);breaks.push_back(-fid);
        const Real selected=selection.restrict_bound(p,fid);
        breaks.push_back(selected);breaks.push_back(-selected);
      }
      const Real direct=integrate([&](Real cosine) {
        auto pair=born_decay(M,cosine);
        pair.first=pair.first.scaled(f1);pair.second=pair.second.scaled(f2);
        if(!cuts.accepts(pair.first,pair.second)) return Real{0};
        const auto observed=[&](const FourVector& v) {
          return selection.cut==SingleLeptonCut::pt_above ? v.pt()>selection.threshold
              : std::abs(v.eta())<selection.threshold;
        };
        return (observed(pair.first)+observed(pair.second))/4.0L;
      },-1,1,{},breaks).value;
      close(isotropic_acceptance(M,f1,f2,cuts,selection),direct,1e-10L,
            "analytic Born cumulative vs literal four-vector selection");
    }
  QEDXSpaceGrid grid;grid.bins=256;
  CoupledDecayOptions options;options.ff_panels=1;
  const std::vector<Real> pt_edges{25,35,M/2,50};
  const std::vector<Real> eta_edges{-2.4L,-.8L,-.4L,0,.4L,.8L,2.4L};
  const CoupledLeptonSpectrum pt(M,cuts,LeptonObservable::transverse_momentum,pt_edges,grid,options);
  const CoupledLeptonSpectrum eta(M,cuts,LeptonObservable::pseudorapidity,eta_edges,grid,options);
  const CoupledLeptonSpectrum asymmetric_eta(M,cuts,LeptonObservable::pseudorapidity,
      {-2.4L,-.4L,.8L,2.4L},grid,options);
  MomentEvolutionParameters p;p.alpha_msbar_at_mu0=0;
  const auto e0=QEDXSpaceEvolution(p,grid);
  const auto rp=pt.match(e0),re=eta.match(e0);
  const Real born=isotropic_acceptance(M,1,1,cuts);
  const Real b35=std::sqrt(1-std::pow(70/M,2));
  close(rp.bins[0].integral.matched,2*(born-b35),1e-13L,"Born pT integral below 35 GeV");
  close(rp.bins[1].integral.matched,2*b35,1e-13L,"Born Jacobian endpoint integral");
  close(rp.bins[2].integral.matched,0,0,"no single muons above M/2");
  close(re.bins[1].integral.matched,std::tanh(.8L)-std::tanh(.4L),1e-13L,"signed Born eta bin");
  close(asymmetric_eta.match(e0).bins[1].integral.matched,std::tanh(.4L)+std::tanh(.8L),
        1e-13L,"asymmetric signed eta bin crossing zero");
  for(const auto& result:{rp,re}) close(result.sum.matched,2*born,1e-13L,"two histogram entries per pair");
  p={};
  for(const auto* spectrum:{&pt,&eta}) {
    const auto central=spectrum->match(QEDXSpaceEvolution(p,grid));
    close(central.matched_closure,0,1e-13L,"single-muon matched closure");
    close(central.fixed_order_closure,0,1e-13L,"single-muon Oalpha closure");
    close(central.coefficient2_closure,0,1e-10L,"single-muon second-coefficient closure");
    auto q=p;q.mu_frag=2*M;
    const auto varied=spectrum->match(QEDXSpaceEvolution(q,grid));
    q.muon_mass/=2;q.mu_0=q.muon_mass;
    const auto low_varied=spectrum->match(QEDXSpaceEvolution(q,grid));
    q.mu_frag=M;
    const auto low_central=spectrum->match(QEDXSpaceEvolution(q,grid));
    for(std::size_t j=0;j<central.bins.size();++j) {
      close(varied.bins[j].integral.continuum_coefficient1,central.bins[j].integral.continuum_coefficient1,
            2e-6L,"single-muon first-order scale cancellation");
      close(low_varied.bins[j].integral.coefficient2-low_central.bins[j].integral.coefficient2,
            varied.bins[j].integral.coefficient2-central.bins[j].integral.coefficient2,3e-6L,
            "single-muon scale shift has no retained alpha2 mass log");
    }
  }
  const auto symmetric=eta.match(QEDXSpaceEvolution(p,grid));
  for(std::size_t j=0;j<symmetric.bins.size();++j)
    close(symmetric.bins[j].integral.matched,symmetric.bins[symmetric.bins.size()-1-j].integral.matched,
          1e-14L,"signed eta reflection in unpolarised rest-frame model");
  // Independent massive real/virtual current establishes the selected leading-
  // power Oalpha target, rather than relying only on histogram closure.
  HardDecayOptions h;h.soft_cutoff=1e-7L;h.collinear_cutoff=1e-8L;
  h.integration={5e-7L,5e-7L,8192};
  for(Real threshold:{40,45}) {
    const SingleLeptonSelection selection{SingleLeptonCut::pt_above,threshold};
    const auto fixed=fiducial_vector_order1(M,cuts,QEDParameters{},h,{},selection);
    close(fixed.geometry_coefficient,vector_muon_pair_geometry(M,cuts,1,1,h,selection),
          3e-6L,"higher-energy-only endpoint contact in fixed and dressed geometry");
  }
  for(SingleLeptonSelection selection:{pt35,{SingleLeptonCut::absolute_eta_below,.5L}}) {
    QEDParameters q;
    const auto lp=fiducial_vector_order1(M,cuts,q,h,{},selection);
    const auto massive=massive_fiducial_vector_order1(M,cuts,q,h,selection);
    close(lp.rate,massive.rate,5e-6L,"selected cumulative massive/leading-power comparison");
  }
  throws<std::invalid_argument>([&]{CoupledLeptonSpectrum(M,cuts,LeptonObservable::transverse_momentum,{25,40},grid);},
      "incomplete pT support");
  throws<std::invalid_argument>([&]{CoupledLeptonSpectrum(M,cuts,LeptonObservable::pseudorapidity,{-2,.5L,2},grid);},
      "incomplete eta support");
  std::cout<<"Single-muon spectra, both-entry normalisation, selected massive and scale checks passed\n";
}
