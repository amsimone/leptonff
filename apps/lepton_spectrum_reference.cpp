#include "leptonff/lepton_spectrum.hpp"
#include "lepton_bins.hpp"
#include <iomanip>
#include <iostream>
#include <limits>

namespace {
using namespace leptonff;
void array(const std::vector<Real>& a) {
  std::cout<<'[';bool first=true;
  for(Real x:a){if(!first)std::cout<<',';first=false;std::cout<<x;}std::cout<<']';
}
void matrix(const std::vector<std::vector<Real>>& a) {
  std::cout<<'[';bool first=true;
  for(const auto& row:a){if(!first)std::cout<<',';first=false;array(row);}std::cout<<']';
}
void rates(const CoupledDecayMatch& r) {
  std::cout<<"{\"born\":"<<r.born<<",\"ff\":"<<r.ff
    <<",\"charge_flipped_pairs\":"<<r.charge_flipped_pairs
    <<",\"muon_projected\":"<<r.muon_projected<<",\"muon_geometry\":"<<r.muon_geometry
    <<",\"photon_projected\":"<<r.photon_projected<<",\"photon_geometry\":"<<r.photon_geometry
    <<",\"muon_hard\":"<<r.a_hard*(r.muon_projected+r.muon_geometry)
    <<",\"photon_hard\":"<<r.a_hard*(r.photon_projected+r.photon_geometry)
    <<",\"overlap_correction\":"<<r.overlap_correction
    <<",\"lattice_coefficient1\":"<<r.lattice_coefficient1
    <<",\"continuum_coefficient1\":"<<r.continuum_coefficient1
    <<",\"coefficient2\":"<<r.coefficient2<<",\"raw\":"<<r.raw
    <<",\"matched\":"<<r.matched<<",\"fixed_order\":"<<r.fixed_order<<'}';
}
void histogram(const CoupledMassSpectrumResult& r,bool pt) {
  std::cout<<"{\"sum\":";rates(r.sum);
  std::cout<<",\"matched_closure\":"<<r.matched_closure
    <<",\"fixed_order_closure\":"<<r.fixed_order_closure
    <<",\"coefficient1_closure\":"<<r.coefficient1_closure
    <<",\"coefficient2_closure\":"<<r.coefficient2_closure<<",\"spectrum\":[";
  bool first=true;
  for(const auto& bin:r.bins) {
    if(!first)std::cout<<',';first=false;
    std::cout<<"{\"lower\":"<<bin.lower<<",\"upper\":"<<bin.upper
      <<(pt?",\"matched_density_per_GeV\":":",\"matched_density_per_eta\":")
      <<bin.integral.matched/(bin.upper-bin.lower)<<",\"integral\":";
    rates(bin.integral);std::cout<<'}';
  }
  std::cout<<"]}";
}
void audit(const CoupledLeptonSpectrum& spectrum) {
  std::cout<<'[';bool first=true;
  for(const auto& op:spectrum.cumulative_operators()) {
    if(!first)std::cout<<',';first=false;
    const auto selection=op.lepton_selection();
    std::cout<<"{\"cut\":"<<static_cast<int>(selection.cut)<<",\"threshold\":"<<selection.threshold
      <<",\"fractions\":";array(op.fractions());
    std::cout<<",\"panel_edges\":";array(op.ff_panel_edges());
    std::cout<<",\"geometry_muon\":";matrix(op.muon_geometry());
    std::cout<<",\"geometry_photon\":";matrix(op.photon_geometry());
    std::cout<<",\"endpoint_muon\":";array(op.muon_endpoint_geometry());
    std::cout<<",\"endpoint_photon_first\":";array(op.photon_first_endpoint_geometry());
    std::cout<<",\"endpoint_photon_second\":";array(op.photon_second_endpoint_geometry());
    std::cout<<",\"geometry_samples\":[";
    // Endpoint and two interior nonidentical fractions test both hard channels.
    bool comma=false;
    for(const auto ij:{std::pair<std::size_t,std::size_t>{0,0},{0,1},{1,2}}) {
      if(ij.second>=op.fractions().size()) continue;
      if(comma)std::cout<<',';comma=true;
      std::cout<<"{\"i\":"<<ij.first<<",\"j\":"<<ij.second
        <<",\"muon\":"<<op.muon_geometry()[ij.first][ij.second]
        <<",\"photon\":"<<op.photon_geometry()[ij.first][ij.second]<<'}';
    }
    std::cout<<"]}";
  }
  std::cout<<']';
}
void point(const char* label,const CoupledLeptonSpectrum& pt,const CoupledLeptonSpectrum& eta,
    const QEDXSpaceEvolution& e,const ReferenceInput& input,const CoupledDecayOptions& options,bool dump) {
  const auto p=pt.match(e,input.integration),y=eta.match(e,input.integration);
  if(std::abs(p.total.matched-y.total.matched)>1e-12L)
    throw std::runtime_error("Single-lepton spectra have inconsistent pair totals");
  std::cout<<"{\"label\":"<<std::quoted(label)<<",\"bins\":"<<e.grid().bins
    <<",\"log_span\":"<<e.grid().log_span<<",\"ff_panels\":"<<options.ff_panels
    <<",\"active_ff_nodes\":"<<pt.pair_operator().fractions().size()
    <<",\"soft_cutoff\":"<<options.geometry.soft_cutoff
    <<",\"collinear_cutoff\":"<<options.geometry.collinear_cutoff
    <<",\"geometry_abs_tol\":"<<options.geometry.integration.absolute_tolerance
    <<",\"geometry_rel_tol\":"<<options.geometry.integration.relative_tolerance
    <<",\"mu_frag\":"<<e.parameters().mu_frag<<",\"mu_0\":"<<e.parameters().mu_0
    <<",\"alpha_msbar_at_mu0\":"<<e.parameters().alpha_msbar_at_mu0
    <<",\"a_initial\":"<<p.total.a_initial<<",\"a_hard\":"<<p.total.a_hard
    <<",\"total\":";rates(p.total);
  std::cout<<",\"weighted_wrap_bound\":"<<e.diagnostics().weighted_wrap_bound
    <<",\"weighted_tail_bound\":"<<e.diagnostics().weighted_tail_bound
    <<",\"compression_moment_residual\":"<<p.total.compression_moment_residual;
  std::cout<<",\"pt\":";histogram(p,true);
  std::cout<<",\"eta\":";histogram(y,false);
  if(dump) {
    const auto& hard=pt.pair_operator();
    const auto n=hard.muon_projection().size();
    std::cout<<",\"hard_muon\":";array(hard.muon_projection());
    std::cout<<",\"hard_photon\":";array(hard.photon_projection());
    for(const auto& entry:{std::make_pair("same",&e.weights()[0][0]),
        std::make_pair("opposite",&e.weights()[0][1]),std::make_pair("gamma",&e.weights()[2][0])}) {
      std::cout<<','<<std::quoted(entry.first)<<':';
      array(std::vector<Real>(entry.second->begin(),entry.second->begin()+n));
    }
    std::cout<<",\"pt_cumulatives\":";audit(pt);
    std::cout<<",\"eta_cumulatives\":";audit(eta);
  }
  std::cout<<'}'<<std::flush;
}
}

int main(int argc,char** argv) {
  try {
    if(argc<3 || argc>6)
      throw std::invalid_argument("Usage: leptonff_lepton_spectrum config lepton_bins [bins [ff_panels [--scan|--dump]]]");
    ReferenceInput input;read_reference_config(argv[1],input);
    const auto edges=read_lepton_bins(argv[2]);
    QEDXSpaceGrid grid;
    CoupledDecayOptions options;
    const auto integer=[](const char* s){std::size_t used=0;auto v=std::stoull(s,&used);
      if(used!=std::string(s).size())throw std::invalid_argument("Invalid integer option");return v;};
    if(argc>=4)grid.bins=integer(argv[3]);
    if(argc>=5)options.ff_panels=integer(argv[4]);
    const std::string mode=argc==6?argv[5]:"";
    if(mode!="" && mode!="--scan" && mode!="--dump")throw std::invalid_argument("Unknown spectrum mode");
    if(mode=="--scan" && (grid.bins>32768 || options.ff_panels>32))
      throw std::invalid_argument("Spectrum scan needs bins<=32768 and panels<=32");
    grid.validate();options.validate();
    options=lepton_spectrum_options(input.mass,input.cuts,edges.pt,edges.eta,options);
    MomentEvolutionParameters p;
    p.alpha_msbar_at_mu0=input.qed.alpha;p.muon_mass=input.qed.muon_mass;
    p.mu_0=input.qed.mu_0;p.mu_frag=input.qed.mu_frag;
    const auto anchor=p;
    std::cout<<std::setprecision(18)<<"{\"version\":\""<<LEPTONFF_VERSION
      <<"\",\"compiler\":"<<std::quoted(__VERSION__)<<",\"real_digits\":"<<std::numeric_limits<Real>::digits
      <<",\"observable\":\"both_muons_in_selected_opposite_charge_pairs\",\"units\":\"dimensionless_integrals; pT_density_per_GeV; signed_eta_density_per_unit\""
      <<",\"bin_convention\":\"two entries per pair; conservative mass cells; pT endpoint M/2; signed eta by rest-frame reflection\""
      <<",\"model\":\"one_family_unpolarised_vector_at_rest_inclusive_distinct_parent_pairs\""
      <<",\"accuracy\":\"Oalpha_hard_matched_coupled_NLL_reference\",\"scheme\":\"MSbar\""
      <<",\"geometry_quadrature\":\"direct_lattice_endpoint_strips_and_contact_aligned_cubic_interior\""
      <<",\"coupling_input\":\"config.alpha is MSbar alpha(mu_0), no OS conversion\""
      <<",\"alpha_msbar_at_mu0\":"<<p.alpha_msbar_at_mu0<<",\"muon_mass\":"<<p.muon_mass
      <<",\"mu_0\":"<<p.mu_0<<",\"mu_frag\":"<<p.mu_frag<<",\"parent_mass\":"<<input.mass
      <<",\"pt_min\":"<<input.cuts.pt_min<<",\"eta_max\":"<<input.cuts.eta_max
      <<",\"mass_min\":"<<input.cuts.mass_min<<",\"mass_max\":"<<input.cuts.mass_max
      <<",\"pt_edges\":";array(edges.pt);std::cout<<",\"eta_edges\":";array(edges.eta);
    std::cout<<",\"ff_fraction_breaks\":";array(options.ff_fraction_breaks);std::cout<<",\"points\":[";
    const QEDXSpaceEvolution e(p,grid);
    const auto make=[&](LeptonObservable observable) {
      return CoupledLeptonSpectrum(input.mass,input.cuts,observable,
          observable==LeptonObservable::transverse_momentum?edges.pt:edges.eta,grid,options);
    };
    auto pt=make(LeptonObservable::transverse_momentum),eta=make(LeptonObservable::pseudorapidity);
    point("baseline",pt,eta,e,input,options,mode=="--dump");
    if(mode=="--scan") {
      options.ff_panels*=2;
      pt=make(LeptonObservable::transverse_momentum);eta=make(LeptonObservable::pseudorapidity);
      std::cout<<',';point("ff_panels_x2",pt,eta,e,input,options,false);
      grid.bins*=2;
      const QEDXSpaceEvolution fine(p,grid);
      pt=make(LeptonObservable::transverse_momentum);eta=make(LeptonObservable::pseudorapidity);
      std::cout<<',';point("grid_x2",pt,eta,fine,input,options,false);
      grid.bins*=2;
      const QEDXSpaceEvolution finest(p,grid);
      pt=make(LeptonObservable::transverse_momentum);eta=make(LeptonObservable::pseudorapidity);
      std::cout<<',';point("grid_x4",pt,eta,finest,input,options,false);
      grid.bins*=2;grid.log_span*=2;
      const QEDXSpaceEvolution span(p,grid);
      pt=make(LeptonObservable::transverse_momentum);eta=make(LeptonObservable::pseudorapidity);
      std::cout<<',';point("span_x2",pt,eta,span,input,options,false);
      options.geometry.soft_cutoff/=4;options.geometry.collinear_cutoff/=4;
      options.geometry.integration.absolute_tolerance/=4;
      options.geometry.integration.relative_tolerance/=4;
      pt=make(LeptonObservable::transverse_momentum);eta=make(LeptonObservable::pseudorapidity);
      std::cout<<',';point("hard_slivers",pt,eta,span,input,options,false);
      for(Real scale:{.5L,2.0L}) {
        p=anchor;p.mu_frag=scale*anchor.mu_frag;
        std::cout<<',';point(scale<1?"scale_half":"scale_double",pt,eta,QEDXSpaceEvolution(p,grid),input,options,false);
      }
      for(Real scale:{.5L,2.0L}) {
        p=anchor;
        auto running=anchor;running.mu_frag=scale*anchor.mu_0;
        p.alpha_msbar_at_mu0=2*std::acos(-1.0L)*evolve_qed_moment(2,running).a_final;
        p.mu_0=running.mu_frag;
        std::cout<<',';point(scale<1?"start_half":"start_double",pt,eta,QEDXSpaceEvolution(p,grid),input,options,false);
      }
    }
    std::cout<<"]}\n";
  }catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}
}
