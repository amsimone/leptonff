// Exercise the actual differential callback without an MCFM integration.
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <sstream>
extern "C" int leptonff_mcfm_evaluate_all(int,const double*,const double*,const double*,const double*,double*,double*);
int main(int argc,char** argv) {
  const int process=argc>1?std::atoi(argv[1]):1;
  if(process!=1 && process!=6 && process!=31)return 2;
  std::ostringstream out;
  out<<std::setprecision(17)<<"[";
  for(int point=0;point<4;++point) {
    const double c=-.7+.45*point,e=(process==31?91.1876+point*10:80.385)/2;
    const double mu[4]={e*std::sqrt(1-c*c),0,e*c,e},other[4]={-mu[0],0,-mu[2],e};
    const auto p3=process==1?other:mu,p4=process==1?mu:other;
    double r[33],s[136],repeat[33],again[136];
    if(leptonff_mcfm_evaluate_all(process,p3,p4,p3,p4,r,s) ||
       leptonff_mcfm_evaluate_all(process,p3,p4,p3,p4,repeat,again))
      throw std::runtime_error("Photonic callback failed");
    for(int k=0;k<33;++k)if(!std::isfinite(r[k]) || r[k]!=repeat[k])
      throw std::runtime_error("Rate finiteness/repeatability failed");
    for(int k=0;k<136;++k)if(!std::isfinite(s[k]) || s[k]!=again[k])
      throw std::runtime_error("Spectrum finiteness/repeatability failed");
    const int widths[3]={process==31?8:6,process==31?5:7,4};
    const int indices[8]={0,15,23,6,25,7,8,24};
    int start=0;
    for(int width:widths) {
      for(int j=0;j<8;++j) {
        double sum=0;for(int k=start;k<start+width;++k)sum+=s[8*k+j];
        if(std::abs(sum-r[indices[j]])>1e-9*(1+std::abs(sum)))
          throw std::runtime_error("Differential/total closure failed");
      }
      start+=width;
    }
    if(point)out<<',';
    out<<"{\"rates\":[";
    for(int k=0;k<33;++k)out<<(k?",":"")<<r[k];
    out<<"],\"spectra\":[";
    for(int k=0;k<136;++k)out<<(k?",":"")<<s[k];
    out<<"]}";
  }
  out<<"]\n";
  std::cout<<"LEPTONFF_JSON="<<out.str();
}
