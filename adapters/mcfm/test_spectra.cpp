// Integration-level checks with external PHOTOS and the actual C/Fortran ABI.
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
extern "C" int leptonff_mcfm_evaluate_all(int,const double*,const double*,const double*,const double*,double*,double*);
int main(int argc,char** argv) {
  const int process=argc>1?std::atoi(argv[1]):1;
  const int count=argc>2?std::atoi(argv[2]):60;
  for(int i=0;i<count;++i) {
    const double c=-.99+1.98*(i+.5)/count,phi=.124*i,e=(process==31?91.1876:80.385)/2;
    const double pt=e*std::sqrt(1-c*c);
    const double mu[4]={pt*std::cos(phi),pt*std::sin(phi),e*c,e};
    const double nu[4]={-mu[0],-mu[1],-mu[2],e};
    const double* p3=process==1?nu:mu;const double* p4=process==1?mu:nu;
    double r[33],again[33],s[136],t[136];
    if(leptonff_mcfm_evaluate_all(process,p3,p4,p3,p4,r,s) ||
       leptonff_mcfm_evaluate_all(process,p3,p4,p3,p4,again,t))throw std::runtime_error("Differential ABI failed");
    for(int k=0;k<136;++k)if(!std::isfinite(s[k]) || s[k]!=t[k])throw std::runtime_error("Spectrum repeatability failed");
    const bool photonic=std::getenv("LEPTONFF_PHOTONIC_LL") && std::string(std::getenv("LEPTONFF_PHOTONIC_LL"))=="1";
    const std::array<int,8> index=photonic?std::array<int,8>{0,15,23,6,25,7,8,24}:std::array<int,8>{0,15,29,30,31,7,8,32};
    const int edges[4]={0,process==31?8:6,13,17};
    for(int obs=0;obs<3;++obs)for(int j=0;j<8;++j) {
      double sum=0;for(int k=edges[obs];k<edges[obs+1];++k)sum+=s[8*k+j];
      if(std::abs(sum-r[index[j]])>2e-9*(1+std::abs(r[index[j]])))throw std::runtime_error("Bin/rate closure failed");
    }
    for(int k=0;k<17;++k) {
      if(std::abs(s[8*k+4]-(s[8*k+2]-s[8*k+1]))>1e-10)throw std::runtime_error("Bin overlap closure failed");
      if(std::abs(s[8*k+7]-(s[8*k+6]-s[8*k+2]))>1e-10)throw std::runtime_error("Paired PHOTOS difference failed");
    }
  }
  std::cout<<"Observed spectra process "<<process<<": "<<count<<" points, all bins finite/repeatable, all eight columns close to each rate\n";
}
