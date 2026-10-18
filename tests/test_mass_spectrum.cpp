#include "leptonff/mass_spectrum.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const Real M=91.1876L;
  FiducialCuts cuts;
  QEDXSpaceGrid grid; grid.bins=256;
  CoupledDecayOptions options; options.ff_panels=1;
  // Include the Born mass as an interior edge: the lower bin must own its
  // atom together with its arbitrarily soft neighbourhood.
  const std::vector<Real> edges{60,80,M,120};
  const CoupledMassSpectrum spectrum(M,cuts,edges,grid,options);
  MomentEvolutionParameters p;
  auto zero=p; zero.alpha_msbar_at_mu0=0;
  const auto identity=spectrum.match(QEDXSpaceEvolution(zero,grid));
  close(identity.bins[0].integral.matched,0,0,"Born outside lower real-tail bin");
  close(identity.bins[1].integral.matched,identity.total.born,1e-14L,
        "Born atom is in bin below an exact M edge");
  close(identity.bins[2].integral.matched,0,0,"no above-endpoint pairs");
  for(Real scale:{.5L,1.0L,2.0L}) {
    p.mu_frag=scale*M;
    const auto e=QEDXSpaceEvolution(p,grid);
    const auto r=spectrum.match(e);
    close(r.matched_closure,0,1e-13L,"matched mass-bin sum closes");
    close(r.fixed_order_closure,0,1e-13L,"fixed-order bin sum closes");
    close(r.coefficient1_closure,0,1e-10L,"first coefficient partitions");
    close(r.coefficient2_closure,0,1e-10L,"second coefficient partitions");
    auto tail=cuts; tail.mass_max=80;
    const QEDParameters q{p.alpha_msbar_at_mu0,p.muon_mass,p.mu_0,p.mu_frag};
    const auto direct=fiducial_vector_order1(M,tail,q,options.geometry);
    close(r.bins[0].integral.continuum_coefficient1,direct.full_coefficient,2e-6L,
          "cumulative bin vs independent direct real-tail Oalpha action");
    close(r.bins[0].integral.fixed_order,direct.rate,3e-9L,"direct real-tail fixed order");
    for(const auto& bin:r.bins) {
      close(bin.integral.a_hard,r.total.a_hard,0,"hard coupling is metadata, not a bin difference");
      close(bin.integral.a_initial,r.total.a_initial,0,"input coupling is shared metadata");
    }
  }
  // The cell measurement must use the same generator/hard subtraction in
  // every bin: no retained alpha^2 mass log may survive a mu_frag variation.
  p={};
  const auto central=spectrum.match(QEDXSpaceEvolution(p,grid));
  p.mu_frag=2*M;
  const auto varied=spectrum.match(QEDXSpaceEvolution(p,grid));
  p.muon_mass/=2;p.mu_0=p.muon_mass;
  const auto low_varied=spectrum.match(QEDXSpaceEvolution(p,grid));
  p.mu_frag=M;
  const auto low_central=spectrum.match(QEDXSpaceEvolution(p,grid));
  for(std::size_t j=0;j<central.bins.size();++j)
    close(low_varied.bins[j].integral.coefficient2-low_central.bins[j].integral.coefficient2,
          varied.bins[j].integral.coefficient2-central.bins[j].integral.coefficient2,2e-6L,
          "cell-integrated mass-bin scale shift contains no retained alpha2 mass log");
  // Restricting a hard operator must preserve the FF quadrature; independently
  // retuning it for each cumulative would spoil a common linear measurement.
  const CoupledFiducialDecay base(M,cuts,grid,options);
  const auto restricted=base.restricted_mass_window(80,120);
  if(base.fractions()!=restricted.fractions() || base.muon_projection()!=restricted.muon_projection())
    throw std::runtime_error("Mass restriction changed shared FF/hard projection");
  throws<std::invalid_argument>([&]{ CoupledMassSpectrum(M,cuts,{60,80,80,120},grid); },"duplicate mass edge");
  throws<std::invalid_argument>([&]{ CoupledMassSpectrum(M,cuts,{60,80,100},grid); },"incomplete mass partition");
  throws<std::invalid_argument>([&]{ base.restricted_mass_window(50,120); },"mass restriction outside support");
  throws<std::invalid_argument>([&]{ CoupledMassSpectrum(M,{25,2.4L,60,M},{60,80,M},grid); },
      "strict fiducial upper cut on Born endpoint remains unsafe");
  std::cout<<"Matched mass spectrum, endpoint ownership and direct tail checks passed\n";
}
