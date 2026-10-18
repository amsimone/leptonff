#include "leptonff/lepton_spectrum.hpp"
#include "leptonff/hard_dressing.hpp"
#include "leptonff/photon_pair.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  const Real M=91.1876L;
  const FiducialCuts cuts;
  QEDXSpaceGrid grid;grid.bins=256;
  CoupledDecayOptions options;options.ff_panels=2;
  options=lepton_spectrum_options(M,cuts,{25,45,50},{-2.4L,0,2.4L},options);
  const Real boundary=std::exp(-grid.spacing());
  options.ff_fraction_breaks.push_back(boundary);
  const SingleLeptonSelection selected{SingleLeptonCut::pt_above,45};
  const auto hard=CoupledFiducialDecay(M,cuts,grid,options).with_lepton_selection(selected);
  MomentEvolutionParameters parameters;
  const QEDXSpaceEvolution e(parameters,grid);
  const auto result=hard.match(e);
  const auto n=hard.muon_projection().size();
  std::vector<Real> at_boundary(n,0);at_boundary[1]=1;
  const auto projected=hard.compress(at_boundary);
  for(std::size_t i=0;i<projected.size();++i) if(hard.fractions()[i]>boundary)
    close(projected[i],0,0,"exact contact belongs to the lower fraction panel");
  // At this spacing no resolved lattice fraction can produce a 45 GeV
  // muon. Resolved/resolved geometry is EXACTLY empty. A panel straddling
  // 90/M could otherwise generate spurious support after compression.
  for(std::size_t i=1;i<hard.fractions().size();++i) {
    if(hard.fractions()[i]>90/M) throw std::runtime_error("Inactive endpoint panel was retained");
    for(std::size_t j=1;j<hard.fractions().size();++j) {
      close(hard.muon_geometry()[i][j],0,0,"muon resolved geometry kinematic support");
      close(hard.photon_geometry()[i][j],0,0,"photon resolved geometry kinematic support");
    }
  }
  Real muon=0,photon=0;
  for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<n;++j) {
    const Real f=e.z(i),g=e.z(j);
    const Real s=e.weights()[0][0][i],o=e.weights()[0][1][i];
    const Real sj=e.weights()[0][0][j],oj=e.weights()[0][1][j];
    muon+=(s*sj+o*oj)*vector_muon_pair_geometry(M,cuts,f,g,options.geometry,selected)/(f*g);
    photon+=2*e.weights()[2][0][i]*(sj+oj)
        *fiducial_vector_photon_pair_hard(M,M,cuts,f,g,options.geometry,selected).geometry/(f*g);
  }
  close(result.muon_geometry,muon,2e-12L,"full signed muon geometry vs literal lattice tensor");
  close(result.photon_geometry,photon,2e-13L,"both photon orientations vs literal lattice tensor");
  close(result.compression_moment_residual,0,2e-12L,"nonuniform FF panels retain cubic moments");
  options.ff_panels=4;
  const auto fine=CoupledFiducialDecay(M,cuts,grid,options).with_lepton_selection(selected);
  const auto refined=fine.match(e);
  close(refined.matched,result.matched,2e-12L,"no endpoint panel dependence on an exactly resolved coarse fixture");
  close(refined.coefficient2,result.coefficient2,2e-9L,"hard alpha2 coefficient uses direct endpoint strips");
  Real previous=0;
  for(Real alpha:{.00025L,.0005L}) {
    parameters.alpha_msbar_at_mu0=alpha;
    const auto r=hard.match(QEDXSpaceEvolution(parameters,grid));
    const Real a=r.a_initial;
    const Real remainder=r.matched-r.born-a*r.continuum_coefficient1-a*a*r.coefficient2;
    if(previous!=0) close(remainder/previous,8,.15L,"corrected geometry has cubic expansion remainder");
    previous=remainder;
  }
  auto invalid=options;invalid.ff_fraction_breaks={1};
  throws<std::invalid_argument>([&]{invalid.validate();},"unphysical FF contact");
  std::cout<<"Direct endpoint strips, threshold-aligned support, signed tensors and expansion passed\n";
}
