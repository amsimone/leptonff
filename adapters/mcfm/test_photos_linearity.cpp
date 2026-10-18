// Statistical check of the signed first-order proposal estimator. Link only
// photos_comparison.o and PHOTOS, run in fresh processes for each scale.
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <stdexcept>
void leptonff_photos_comparison(int,const double*,const double*,const double*,const double*,double*);
int main(int argc,char** argv) {
  if(argc!=4)return 2;
  const int process=std::atoi(argv[1]),count=std::atoi(argv[2]);
  std::ofstream out(argv[3]);
  for(int i=0;i<count;++i) {
    const double c=-.999+1.998*(i+.5)/count,phi=.124*i,e=(process==31?91.1876:80.385)/2;
    const double pt=e*std::sqrt(1-c*c);
    const double mu[4]={pt*std::cos(phi),pt*std::sin(phi),e*c,e},nu[4]={-mu[0],-mu[1],-mu[2],e};
    double r[4];
    leptonff_photos_comparison(process,process==1?nu:mu,process==1?mu:nu,mu,nu,r);
    if(!std::isfinite(r[2]))throw std::runtime_error("Nonfinite proposal correction");
    out<<r[2]<<'\n';
  }
}
