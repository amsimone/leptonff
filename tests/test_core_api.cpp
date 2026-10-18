#include "leptonff/core.hpp"
#include "leptonff/fragmentation_c.h"
#include "check.hpp"
#include <cstring>

double coefficient(double x,void*) { return x*(1-x); }
int main() {
  using namespace leptonff;
  if(std::strcmp(leptonff_version(),LEPTONFF_EXPECTED_VERSION) || std::strcmp(version(),leptonff_version()))
    throw std::runtime_error("Native version reporting mismatch");
  leptonff_fragmentation_config c;
  leptonff_fragmentation_config_default(&c);
  leptonff_fragmentation_context* context=nullptr;
  if(leptonff_fragmentation_create(&c,&context)) throw std::runtime_error(leptonff_last_error());
  double value=0;
  const auto ff=fragmentation_order1(Species::muon_minus,Species::muon_minus,{});
  for(double N : {1.,2.,3.}) {
    if(leptonff_fragmentation_moment(context,0,0,N,&value)) throw std::runtime_error(leptonff_last_error());
    close(value,ff.moment(N),1e-12L,"C/C++ fragmentation moments");
  }
  if(leptonff_fragmentation_integral(context,0,0,.9,1,&value)) throw std::runtime_error(leptonff_last_error());
  close(value,ff.integral(.9L,1),1e-10L,"C/C++ endpoint bin");
  if(leptonff_fragmentation_convolve(context,0,0,.4,coefficient,nullptr,&value)) throw std::runtime_error(leptonff_last_error());
  close(value,ff.convolve([](Real x){return x*(1-x);},.4L),1e-10L,"C callback convolution");
  if(!leptonff_fragmentation_moment(context,99,0,2,&value) || !std::strlen(leptonff_last_error()))
    throw std::runtime_error("C API failed to reject invalid species");
  leptonff_fragmentation_destroy(context);
  c.model=LEPTONFF_FF_COUPLED_NLL;
  if(leptonff_fragmentation_create(&c,&context)) throw std::runtime_error(leptonff_last_error());
  if(leptonff_fragmentation_moment(context,0,0,2,&value)) throw std::runtime_error(leptonff_last_error());
  close(value,evolve_qed_moment(2).fragmentation[0][0],1e-11L,"C/C++ coupled NLL");
  if(!leptonff_fragmentation_integral(context,0,0,.1,1,&value))
    throw std::runtime_error("Moment-only context accepted an x-space bin");
  leptonff_fragmentation_destroy(context);
  c.alpha=0; c.grid_bins=128;
  if(leptonff_fragmentation_create(&c,&context)) throw std::runtime_error(leptonff_last_error());
  if(leptonff_fragmentation_integral(context,0,0,.8,1,&value)) throw std::runtime_error(leptonff_last_error());
  close(value,1,1e-12L,"prepared coupled endpoint identity");
  if(leptonff_fragmentation_integral(context,0,0,.2,.8,&value)) throw std::runtime_error(leptonff_last_error());
  close(value,0,1e-12L,"prepared coupled interior identity");
  if(!leptonff_fragmentation_integral(context,0,0,0,1,&value))
    throw std::runtime_error("Coupled x-space action accepted zero lower support");
  leptonff_fragmentation_destroy(context);
  c.abi_version=999;
  if(!leptonff_fragmentation_create(&c,&context) || context)
    throw std::runtime_error("Invalid ABI accepted or dangling output handle");
}
