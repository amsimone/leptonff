// SPDX-License-Identifier: MIT
// Inclusive mass only: conservative production reconstruction and paired PHOTOS.
#include "mass_radiator.hpp"
#include "photos_decay.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <chrono>
#include <sstream>
using namespace leptonff;
namespace {
struct Source {Real lo,hi,weight,error;std::array<Real,5> polynomial{};
  std::size_t first=0;std::array<std::array<Real,5>,5> basis{};};
std::vector<std::pair<Real,Real>> gauss(unsigned n) {
  const Real pi=std::acos(-1.L);std::vector<std::pair<Real,Real>> out;
  for(unsigned i=0;i<n;++i){Real x=std::cos(pi*(i+.75L)/(n+.5L)),derivative=0;
    for(int k=0;k<100;++k){Real p=1,prev=0;for(unsigned j=1;j<=n;++j){const Real old=p;p=((2*j-1)*x*p-(j-1)*prev)/j;prev=old;}
      derivative=n*(x*p-prev)/(x*x-1);const Real delta=p/derivative;x-=delta;if(std::abs(delta)<1e-18L)break;}
    out.push_back({(1+x)/2,1/((1-x*x)*derivative*derivative)});
  }return out;
}
Real density(const Source& s,Real M) {
  const Real t=(M-(s.lo+s.hi)/2)/(s.hi-s.lo);
  Real v=0;for(int k=4;k>=0;--k)v=v*t+s.polynomial[k];
  if(v<0)throw std::runtime_error("Negative reconstructed production density");return v;
}
}
int main(int argc,char** argv) {
 try {
  if(argc!=8)throw std::runtime_error("Usage: mass_precision production.tsv output.json grid_bins nodes photos_trials seed alpha");
  if(std::ifstream(argv[2]))throw std::runtime_error("Output exists");
  const auto start=std::chrono::steady_clock::now();
  std::ifstream in(argv[1]);std::vector<Source> source;Source row;
  while(in>>row.lo>>row.hi>>row.weight>>row.error)source.push_back(row);
  if(!in.eof() || source.size()<5)throw std::runtime_error("Malformed production table");
  const Real h=source[0].hi-source[0].lo;
  if(source.front().lo!=66 || source.back().hi!=116)throw std::runtime_error("Requires parent masses 66..116");
  Real normalization=0;
  for(std::size_t j=0;j<source.size();++j){auto& s=source[j];
    if(std::abs(s.hi-s.lo-h)>1e-10 || s.weight<=0 || s.error<0 || (j && source[j-1].hi!=s.lo))throw std::runtime_error("Invalid uniform source grid");
    // Quartic density whose integrals reproduce five neighbouring native
    // bins. Reconstruction preserves each input integral exactly.
    const std::size_t first=std::min(source.size()-5,j>1?j-2:0);
    s.first=first;Real matrix[5][10]{};
    for(int k=0;k<5;++k){const Real center=Real(first+k)-j;
      for(int power=0;power<5;++power)matrix[k][power]=(std::pow(center+.5L,power+1)-std::pow(center-.5L,power+1))/(power+1);
      matrix[k][5+k]=1;}
    for(int k=0;k<5;++k){int pivot=k;for(int row=k+1;row<5;++row)if(std::abs(matrix[row][k])>std::abs(matrix[pivot][k]))pivot=row;
      for(int col=k;col<10;++col)std::swap(matrix[k][col],matrix[pivot][col]);
      const Real d=matrix[k][k];for(int col=k;col<10;++col)matrix[k][col]/=d;
      for(int row=0;row<5;++row)if(row!=k){const Real f=matrix[row][k];for(int col=k;col<10;++col)matrix[row][col]-=f*matrix[k][col];}}
    for(int k=0;k<5;++k)for(int j=0;j<5;++j){s.basis[k][j]=matrix[k][5+j]/h;s.polynomial[k]+=s.basis[k][j]*source[first+j].weight;}
    normalization+=s.weight;
  }
  LLGrid grid;grid.bins=std::stoul(argv[3]);grid.log_span=32;grid.soft_log_cutoff=grid.spacing()/2;
  const unsigned nodes=std::stoul(argv[4]);const auto draws=std::stoull(argv[5]),seed=std::stoull(argv[6]);
  if(nodes<2 || nodes>1024)throw std::runtime_error("Quadrature node count must be 2..1024");
  QEDParameters q;q.alpha=std::stold(argv[7]);q.mu_frag=91.1876L;
  leptonff_dyturbo::MassRadiator radiator(q,grid);
  constexpr unsigned bins=265; // underflow + 264 visible bins, 0.25 GeV
  std::array<std::array<Real,3>,bins> mean{};
  std::array<std::array<Real,9>,bins> covariance{};
  std::vector<std::array<std::array<Real,3>,bins>> source_response(source.size());
  std::array<Real,bins> photos_mean{},photos_variance{};
  std::array<Real,bins> control_mean{},control_variance{};
  Real peak_photos_delta=0,peak_photos_variance=0;
  const bool control=std::getenv("LEPTONFF_PHOTOS_CONTROL")!=nullptr;
  const std::string control_kind=control?std::getenv("LEPTONFF_PHOTOS_CONTROL"):"none";
  const auto environment=[](const char* key,Real fallback){const char* value=std::getenv(key);return value?std::stold(value):fallback;};
  const Real photos_cutoff=environment("LEPTONFF_PHOTOS_CUTOFF",1e-5L),photos_precision=environment("LEPTONFF_PHOTOS_EXP_EPS",1e-4L);
  std::ostringstream transfer;transfer<<std::setprecision(18)<<'[';
  const auto rule=gauss(nodes);
  for(const auto& s:source){
    const auto segments=std::max(1u,unsigned(std::llround(h/.25L)));const Real segment_width=h/segments;
    if(std::abs(.25L/segment_width-std::round(.25L/segment_width))>1e-10L)throw std::runtime_error("Source bins must align with observed 0.25 GeV edges");
    for(unsigned segment=0;segment<segments;++segment)for(const auto [t,w]:rule){
      // Smooth endpoint map reduces the integrable soft-endpoint log.
      const Real u=t*t*(3-2*t),jac=6*t*(1-t),M=s.lo+segment_width*(segment+u);
      const Real local=(M-(s.lo+s.hi)/2)/h;
      std::array<Real,5> basis{};
      for(int j=0;j<5;++j){for(int k=4;k>=0;--k)basis[j]=basis[j]*local+s.basis[k][j];basis[j]*=w*jac*segment_width;}
      std::array<Real,3> accounted{};
      for(unsigned b=1;b<bins;++b){const auto r=radiator.bin(M,50+(b-1)*.25L,50+b*.25L);
        for(int a=0;a<3;++a){for(int j=0;j<5;++j)source_response[s.first+j][b][a]+=basis[j]*r[a];accounted[a]+=r[a];}}
      for(int a=0;a<3;++a)for(int j=0;j<5;++j)source_response[s.first+j][0][a]+=basis[j]*(1-accounted[a]);
    }
  }
  for(std::size_t j=0;j<source.size();++j)for(unsigned b=0;b<bins;++b)for(int a=0;a<3;++a){
    mean[b][a]+=source[j].weight*source_response[j][b][a];
    for(int c=0;c<3;++c)covariance[b][3*a+c]+=source[j].error*source[j].error*source_response[j][b][a]*source_response[j][b][c];}
  std::uint64_t actual=0;
  if(draws){
    leptonff_dyturbo::PhotosDecay photos(q.alpha,photos_cutoff,photos_precision);
    std::unique_ptr<leptonff_dyturbo::PhotosDecay> tight;
    if(control)tight=std::make_unique<leptonff_dyturbo::PhotosDecay>(q.alpha,
      control_kind=="exponentiation"?1e-5:1e-7,control_kind=="cutoff"?1e-4:1e-7);
    std::mt19937_64 generator(seed);const auto uniform=[&](){return Real((generator()>>12)+.5L)/4503599627370496.L;};
    const auto bin=[](Real m){return m<50?0u:std::min(bins-1,unsigned((m-50)*4)+1);};
    for(std::size_t j=0;j<source.size();++j){const auto& s=source[j];
      // Stratify in parent mass; mix cross-section and uniform allocation.
      const auto n=std::max(2ULL,static_cast<unsigned long long>(std::ceil(draws*(.8L*s.weight/normalization+.2L/source.size()))));
      std::array<Real,bins> sum{},square{};
      std::array<Real,bins> control_sum{},control_square{};
      Real peak_sum=0,peak_square=0;
      std::array<std::array<Real,5>,bins> moments{};
      std::array<std::array<Real,15>,bins> moment_products{};
      for(std::uint64_t k=0;k<n;++k){const Real M=s.lo+h*uniform(),v=h*density(s,M);
        const BornDecayPoint point{{M/2,M/2,0,0},{M/2,-M/2,0,0},generator()};
        const auto pair=photos.draw(SpinDecay::neutral_vector,point,point.radiation_seed);
        const auto before=bin(M),after=bin(invariant_mass(pair.first,pair.second));
        const auto peak=[](unsigned b){return b>=161 && b<=168;};
        const Real delta=v*(int(peak(after))-int(peak(before)));peak_sum+=delta;peak_square+=delta*delta;
        if(control){const auto refined=tight->draw(SpinDecay::neutral_vector,point,point.radiation_seed);
          const auto changed=bin(invariant_mass(refined.first,refined.second));
          if(changed!=after){control_sum[changed]+=v;control_sum[after]-=v;control_square[changed]+=v*v;control_square[after]+=v*v;}}
        if(before!=after){sum[after]+=v;sum[before]-=v;square[after]+=v*v;square[before]+=v*v;
          std::array<Real,9> powers{};powers[0]=1;const Real t=(M-(s.lo+s.hi)/2)/h;
          for(int a=1;a<9;++a)powers[a]=powers[a-1]*t;
          for(int a=0;a<5;++a){moments[after][a]+=powers[a];moments[before][a]-=powers[a];}
          int index=0;for(int a=0;a<5;++a)for(int c=a;c<5;++c){moment_products[after][index]+=powers[a+c];moment_products[before][index++]+=powers[a+c];}
        }
      }
      for(unsigned b=0;b<bins;++b){photos_mean[b]+=sum[b]/n;photos_variance[b]+=(square[b]-sum[b]*sum[b]/n)/(n*(n-1));
        control_mean[b]+=control_sum[b]/n;control_variance[b]+=(control_square[b]-control_sum[b]*control_sum[b]/n)/(n*(n-1));}
      peak_photos_delta+=peak_sum/n;peak_photos_variance+=(peak_square-peak_sum*peak_sum/n)/(Real(n)*(n-1));
      actual+=n;
      if(j)transfer<<',';transfer<<"{\"lower\":"<<s.lo<<",\"upper\":"<<s.hi<<",\"trials\":"<<n<<",\"polynomial\":[";
      for(int a=0;a<5;++a)transfer<<(a?",":"")<<s.polynomial[a];transfer<<"],\"response\":[";
      for(unsigned b=0;b<bins;++b){if(b)transfer<<',';transfer<<'[';
        for(int a=0;a<5;++a)transfer<<(a?",":"")<<moments[b][a]/n;
        int index=0;for(int a=0;a<5;++a)for(int c=a;c<5;++c)transfer<<','<<(moment_products[b][index++]-moments[b][a]*moments[b][c]/n)/(Real(n)*(n-1));
        transfer<<']';}transfer<<"]}";
    }
    for(unsigned b=0;b<bins;++b)photos_mean[b]+=mean[b][0];
  }
  std::ofstream out(argv[2]);out<<std::setprecision(18)<<"{\"schema\":\"leptonff-mass-precision-1\",\"grid\":"<<grid.bins
    <<",\"nodes\":"<<nodes<<",\"photos_trials\":"<<actual<<",\"seed\":"<<seed<<",\"alpha\":"<<q.alpha<<",\"photos_full_me\":false,\"bins\":[";
  for(unsigned b=0;b<bins;++b){out<<(b?",":"")<<"{\"lower\":"<<(b?50+(b-1)*.25:0)<<",\"upper\":"<<(b?50+b*.25:50)<<",\"mean\":[";
    for(int a=0;a<3;++a)out<<(a?",":"")<<mean[b][a];out<<"],\"covariance\":[";
    for(int a=0;a<9;++a)out<<(a?",":"")<<covariance[b][a];out<<"],\"photos\":"<<photos_mean[b]<<",\"photos_variance\":"<<photos_variance[b]
      <<",\"photos_control_delta\":"<<control_mean[b]<<",\"photos_control_variance\":"<<control_variance[b]<<'}';}
  out<<"],\"transfer_moments\":"<<transfer.str()<<"],\"photos_version\":\""<<leptonff_dyturbo::PhotosDecay::version()<<"\",\"photos_cutoff\":"<<photos_cutoff<<",\"photos_exponentiation_precision\":"<<photos_precision<<",\"photos_control\":"<<(control?"true":"false")
    <<",\"photos_control_kind\":\""<<control_kind<<"\",\"peak_photos_delta\":"<<peak_photos_delta<<",\"peak_photos_variance\":"<<peak_photos_variance
    <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
  if(!out)throw std::runtime_error("Cannot write result");
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
}
