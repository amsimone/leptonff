// Optional independent tree-amplitude check against RECOLA 1.5.1.
#include "recola.h"
#include "leptonff/spin_decay.hpp"
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
using namespace leptonff;
void put(double* p,const FourVector& v) {p[0]=v.energy;p[1]=v.px;p[2]=v.py;p[3]=v.pz;}
int main(int argc,char** argv) {
 if(argc!=2)return 2;
 Recola::set_output_file_rcl("recola.log");
 Recola::set_print_level_squared_amplitude_rcl(0);
 Recola::set_pole_mass_w_rcl(80.385,0);Recola::set_pole_mass_z_rcl(91.1876,0);
 Recola::set_light_muon_rcl();Recola::use_alpha0_scheme_rcl(1/137.035999084);
 const std::string pol[3]={"[-]","[0]","[+]"};
 for(int kind=0;kind<3;++kind) {
  const std::string parent=kind==0?"W-":kind==1?"W+":"Z";
  const std::string daughters=kind==0?"mu- nu_mu~":kind==1?"mu+ nu_mu":"mu- mu+";
  Recola::define_process_rcl(10*kind+1,parent+" -> "+daughters,"LO");
  for(int h=0;h<3;++h)Recola::define_process_rcl(10*kind+h+2,parent+pol[h]+" -> "+daughters+" A","LO");
 }
 Recola::generate_processes_rcl();
 std::ofstream out(argv[1]);out<<std::setprecision(17)<<"kind,helicity,z,c,orientation,twist,recola_coefficient,trace_coefficient,relative_difference\n";
 Real maximum=0;
 for(int kind=0;kind<3;++kind) {
  const auto decay=kind==0?SpinDecay::wminus:kind==1?SpinDecay::wplus:SpinDecay::neutral_vector;
  const Real M=kind==2?91.1876L:80.385L;
  const FourVector parent{std::sqrt(M*M+100),0,0,10};
  double bornp[3][4],amplitude[2];put(bornp[0],parent);
  put(bornp[1],rest_frame_boost({M/2,M/2,0,0},parent,1));
  put(bornp[2],rest_frame_boost({M/2,-M/2,0,0},parent,1));
  Recola::compute_process_rcl(10*kind+1,bornp,"LO",amplitude);
  const Real born=amplitude[0],e2=4*std::acos(-1.L)/137.035999084L;
  for(Real z:{.2L,.6L,.95L})for(Real c:{-.8L,.3L,.9L})for(Real orientation:{-.6L,.4L})for(Real twist:{.1L,1.1L,2.4L}) {
   const auto point=spin_radiative_point(decay,M,0,z,c);
   const auto rotation=decay_rotation(orientation,.7L,twist);
   const auto tensor=rotate_spin_tensor(point.real,rotation);
   double p[4][4];put(p[0],parent);
   put(p[1],rest_frame_boost(rotate_momentum(point.first,rotation),parent,1));
   put(p[2],rest_frame_boost(rotate_momentum(point.second,rotation),parent,1));
   put(p[3],rest_frame_boost(rotate_momentum(point.photon,rotation),parent,1));
   for(int h=0;h<3;++h) {
    SpinTensor rho{};
    if(h==1)rho[2][2]=1;
    else {rho[0][0]=rho[1][1]=.5L;rho[0][1]=std::complex<Real>{0,Real(h-1)/2};rho[1][0]=std::conj(rho[0][1]);}
    Recola::compute_process_rcl(10*kind+h+2,p,"LO",amplitude);
    Real ref=amplitude[0];
    if(kind==2) {
     for(int j=0;j<4;++j)std::swap(p[1][j],p[2][j]);
     Recola::compute_process_rcl(10*kind+h+2,p,"LO",amplitude);ref=(ref+amplitude[0])/2;
     for(int j=0;j<4;++j)std::swap(p[1][j],p[2][j]);
    }
    ref=ref/born*M*M*(1-z)/(4*e2);
    const Real ours=spin_contract(rho,tensor),relative=(ours-ref)/std::max(std::abs(ref),Real{1e-12L});
    maximum=std::max(maximum,std::abs(relative));
    out<<kind<<','<<h-1<<','<<z<<','<<c<<','<<orientation<<','<<twist<<','<<ref<<','<<ours<<','<<relative<<'\n';
   }
  }
 }
 std::cout<<"maximum relative polarised squared-amplitude difference="<<maximum<<'\n';
 Recola::reset_recola_rcl();
 if(maximum>1e-7L)throw std::runtime_error("Independent polarised real-current check failed");
}
