#include "leptonff/spin_decay.hpp"
#include "leptonff/vector_current.hpp"
#include "leptonff/fragmentation.hpp"
#include "check.hpp"
#include <algorithm>

int main() {
  using namespace leptonff;
  const Real M=80.385L,m=.1056583755L,pi=std::acos(-1.L);
  const SpinTensor rho{{{{.6L,{.03L,.04L},{-.02L,.01L}}},
                       {{{.03L,-.04L},.3L,{.02L,-.03L}}},
                       {{{-.02L,-.01L},{.02L,.03L},.1L}}}};
  // Sphere cubature exact through degree 4; tests all nine spin components,
  // including imaginary coherences and both CP-conjugate Born analysers.
  const std::array<Real,3> nodes{-std::sqrt(3.L/5),0,std::sqrt(3.L/5)};
  const std::array<Real,3> weights{5.L/18,4.L/9,5.L/18};
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector}) {
    SpinTensor reconstructed{};
    for(int k=0;k<3;++k)for(int j=0;j<8;++j) {
      const Real c=nodes[k],s=std::sqrt(1-c*c),phi=2*pi*(j+.5L)/8;
      const std::array<Real,3> n{s*std::cos(phi),s*std::sin(phi),c};
      const auto b=spin_born_tensor(kind,n),dual=spin_angular_projector(kind,n);
      const Real w=weights[k]*spin_contract(rho,b)/8;
      for(int a=0;a<3;++a)for(int d=0;d<3;++d)reconstructed[a][d]+=w*dual[a][d];
      close(spin_trace(dual),1,1e-14L,"dual projector trace");
    }
    for(int a=0;a<3;++a)for(int d=0;d<3;++d) {
      close(reconstructed[a][d].real(),rho[a][d].real(),1e-13L,"spin tomography symmetric components");
      close(reconstructed[a][d].imag(),kind==SpinDecay::neutral_vector?0:rho[a][d].imag(),1e-13L,"spin tomography parity components");
    }
    for(Real z:{.2L,.6L,.9L}) {
      const auto p=spin_radiative_point(kind,M,m,z,.3L);
      close(p.first.energy+p.second.energy+p.photon.energy,M,1e-12L,"spin decay energy");
      close(invariant_mass(p.first,p.second),M*std::sqrt(z),1e-12L,"spin decay pair mass");
      close(spin_photon_ward_residual(kind,M,m,z,.3L),0,1e-12L,"radiative photon Ward identity");
      const Real r=m*m/(M*M*z),beta=kind==SpinDecay::neutral_vector?std::sqrt(1-4*r):(1-r)/(1+r);
      const Real edge=std::atanh(beta);
      const Real numeric=integrate([&](Real y) {
        return spin_trace(spin_radiative_point(kind,M,m,z,std::tanh(y)/beta).real)/3
          /(beta*std::pow(std::cosh(y),2));
      },-edge,edge,{1e-6L,1e-7L,8192}).value;
      const Real P=(1+z*z)/(1-z),L=2*std::log(M/m);
      const Real expected=kind==SpinDecay::neutral_vector?2*P*(L+std::log(z)-1)
        :P*(L+std::log(z)-1)-2*(1+z+z*z)/(3*(1-z));
      close(numeric,expected,.0006L,"massive trace vs independent leading-power real spectrum");
      if(kind==SpinDecay::neutral_vector)
        close(spin_trace(p.real)/3,vector_real_density(M,m,z,.3L),1e-12L,"independent vector current trace");
      const auto rotation=decay_rotation(.4L,.7L,.9L);
      close(spin_trace(rotate_spin_tensor(p.real,rotation)),spin_trace(p.real),1e-12L,"tensor rotation trace");
      close(spin_contract(rotate_spin_tensor(rho,rotation),rotate_spin_tensor(p.real,rotation)),
            spin_contract(rho,p.real),1e-12L,"spin contraction covariance");
      const FourVector parent{std::sqrt(M*M+1300),10,20,30};
      const auto lab=rest_frame_boost(p.first,parent,1),back=rest_frame_boost(lab,parent,-1);
      close(back.energy,p.first.energy,1e-12L,"boost inverse energy");
      close(back.pz,p.first.pz,1e-12L,"boost inverse momentum");
    }
    // The soft current and the lepton-collinear current must retain the
    // complete Born spin matrix, including its parity-odd imaginary part.
    for(bool soft:{true,false}) {
      const Real z=soft?1-1e-7L:.6L,c=soft?.3L:1-1e-8L;
      const auto p=spin_radiative_point(kind,M,0,z,c);
      const FourVector cluster{p.first.energy+p.photon.energy,p.first.px,p.first.py,p.first.pz+p.photon.pz};
      const Real norm=std::hypot(cluster.pt(),cluster.pz);
      const auto b=spin_born_tensor(kind,{cluster.px/norm,cluster.py/norm,cluster.pz/norm});
      for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        close(std::abs(p.real[i][j]/(spin_trace(p.real)/3)-b[i][j]),0,soft?1e-6L:1e-3L,"soft/collinear spin factorisation");
    }
  }
  // A mass-window measurement has zero geometric remainder for an unpolarised
  // source at EVERY point. Its complete H1 is the analytic projected kernel.
  SpinTensor unpolarised{};for(int i=0;i<3;++i)unpolarised[i][i]=1.L/3;
  for(Real z:{.3L,.7L,.95L})for(Real c:{-.7L,.2L,.9L}) {
    const auto g=decay_hard_geometry(SpinDecay::neutral_vector,M,z,c,decay_rotation(.2L,.3L,.6L),
      {M,0,0,0},[](const FourVector& a,const FourVector& b){return invariant_mass(a,b)>50?1:0;});
    close(spin_contract(unpolarised,g),0,1e-11L,"inclusive-angle mass-window geometry");
  }
  // Independent inclusive-angle moment and factorisation-scale cancellation.
  for(auto kind:{SpinDecay::wminus,SpinDecay::neutral_vector})for(Real scale:{M/2,M,2*M}) {
    const int legs=kind==SpinDecay::neutral_vector?2:1;
    QEDParameters q;q.mu_frag=scale;
    const Real hard=legs*integrate([&](Real z){
      return projected_decay_hard_kernel(kind,z,M,scale)*(z-1);
    },0,1).value;
    const Real ff=legs*fragmentation_coefficient_order1(Species::muon_minus,Species::muon_minus,q).moment(2);
    const Real L=2*std::log(M/q.muon_mass);
    close(hard+ff,legs==1?11.L/3-4*L/3:44.L/9-8*L/3,2e-8L,
          "unitary physical moment and fragmentation-scale cancellation");
  }
  // Reference subtraction is an orientation-independent zero-integral control
  // variate. Fully inclusive radiation then cancels exactly at every point.
  for(auto kind:{SpinDecay::wminus,SpinDecay::wplus,SpinDecay::neutral_vector}) {
    const auto g=decay_hard_geometry(kind,M,.7L,.2L,decay_rotation(.2L,.3L,.6L),
      {M,0,0,0},[](const FourVector&,const FourVector&){return 1;},1);
    for(auto row:g)for(auto x:row)close(std::abs(x),0,0,"unitarity control variate");
  }
  throws<std::invalid_argument>([&]{spin_born_tensor(SpinDecay::wminus,{1,1,0});},"nonunit spin direction");
  throws<std::invalid_argument>([&]{spin_radiative_point(SpinDecay::wminus,M,0,.5L,1);},"massless collinear endpoint");
  std::cout<<"Polarised W/Z real currents, Ward, spectra, spin tomography and limits passed\n";
}
