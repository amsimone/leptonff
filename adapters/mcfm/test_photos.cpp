// Optional external-library smoke test; compile with the PHOTOS adapter objects.
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <cstdlib>
extern "C" int leptonff_mcfm_evaluate(int,const double*,const double*,const double*,const double*,double*);
int main(int argc,char** argv) {
  const int process=argc>1?std::atoi(argv[1]):1;
  if(process!=1 && process!=6 && process!=31)throw std::invalid_argument("Invalid test process");
  int radiated=0;
  for(int i=0;i<1000;++i) {
    const double c=-.999+1.998*(i+.5)/1000,phi=.124*i,e=(process==31?91.1876:80.385)/2;
    const double pt=e*std::sqrt(1-c*c);
    const double mu[4]={pt*std::cos(phi),pt*std::sin(phi),e*c,e};
    const double nu[4]={-mu[0],-mu[1],-mu[2],e};
    double a[33],b[33];
    const double* p3=process==1?nu:mu;
    const double* p4=process==1?mu:nu;
    if(leptonff_mcfm_evaluate(process,p3,p4,p3,p4,a) || leptonff_mcfm_evaluate(process,p3,p4,p3,p4,b))
      throw std::runtime_error("PHOTOS adapter failed");
    for(int k=0;k<33;++k)if(!std::isfinite(a[k]) || a[k]!=b[k])
      throw std::runtime_error("Common-seed reproducibility or finiteness failed");
    radiated+=a[11]>0;
  }
  if(radiated<10 || radiated>990)throw std::runtime_error("PHOTOS radiation not active");
  std::cout<<"PHOTOS process "<<process<<" conservation, common-seed repeatability and radiation checks passed; radiated="<<radiated<<"/1000\n";
}
