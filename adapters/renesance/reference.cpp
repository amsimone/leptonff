// SPDX-License-Identifier: MIT
// Independent fixed-order reference: upstream SANC Born/SV/real amplitudes.
#include "Info.h"
#include "Settings.h"
#include "Process_pp.h"
#include "ModulesInterface.h"
#include "PartonDidtributionsLHAPDF.h"
#include "defaults.hpp"
#include "photos_decay.hpp"
#include "pythia_decay.hpp"
#include "mass_radiator.hpp"
#include "leptonff/decay_setup.hpp"
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
extern "C" {
void ffini_();
void cc_ha_1413_1615_1spr_(double*,double*,double*);
void cc_ha_1314_1516_1spr_(double*,double*,double*);
void nc_fsr_ha_1313_1616_1spr_(double*,double*,double*);
void nc_fsr_ha_1414_1616_1spr_(double*,double*,double*);
}
namespace SANC {Settings::Settings(){} Settings::~Settings(){}}
namespace {
using namespace leptonff;
using namespace SANCmodulesVar;
class Inputs final:public SANC::Settings {
public:
  std::map<std::string,std::string> v=sanc_defaults();
  void init() override{} void print_diff() override{} void print() override{} void validate() override{}
  bool get_b(std::string k) override{return v.at(k)=="true";}
  long get_l(std::string k) override{return std::stol(v.at(k));}
  double get_d(std::string k) override{return std::stod(v.at(k));}
  std::string get_s(std::string k) override{return v.at(k);}
  void set(std::string k,bool x) override{v[k]=x?"true":"false";}
  void set(std::string k,long x) override{v[k]=std::to_string(x);}
  void set(std::string k,double x) override{std::ostringstream s;s<<std::setprecision(17)<<x;v[k]=s.str();}
  void set(std::string x) override{const auto p=x.find('=');if(p==x.npos)throw std::runtime_error("Expected key=value");v[x.substr(0,p)]=x.substr(p+1);}
};
struct Bin {std::string name;double lo,hi;};
constexpr unsigned channels=9;
using Values=std::array<long double,channels>;
FourVector read(const SANC::Vec4& p){return {p.e(),p.px(),p.py(),p.pz()};}
struct Histograms {
  bool neutral;std::vector<Bin> bins;
  std::vector<Values> sums;
  std::vector<std::array<long double,channels*channels>> products;
  Histograms(bool z):neutral(z) {
    bins.push_back({"inclusive",0,0});bins.push_back({"fiducial",0,0});
    const auto add=[&](std::string s,double a,double b,double h){for(double x=a;x<b-h/2;x+=h)bins.push_back({s,x,x+h});};
    if(z)add("mll",50,116,.25);
    else{add("ptl",25,55,.25);add("absetal",0,2.4,.05);add("mt",50,100,.5);}
    sums.resize(bins.size());products.resize(bins.size());
  }
  std::vector<long double> measure(const FourVector& a,const FourVector& b)const {
    std::vector<long double> f(bins.size());f[0]=1;
    const auto M=invariant_mass(a,b);
    const bool pass=neutral?(M>=50 && M<116):(a.pt()>25 && std::abs(a.eta())<2.4);
    if(!pass)return f;f[1]=1;
    const auto mt=std::sqrt(std::max(Real{0},2*(a.pt()*b.pt()-a.px*b.px-a.py*b.py)));
    for(std::size_t i=2;i<bins.size();++i){const auto& bin=bins[i];
      const Real x=bin.name=="mll"?M:bin.name=="ptl"?a.pt():bin.name=="mt"?mt:std::abs(a.eta());
      f[i]=x>=bin.lo && x<bin.hi;
    }return f;
  }
  void add(const std::vector<Values>& values){
    for(std::size_t i=0;i<bins.size();++i)for(unsigned a=0;a<channels;++a){
      const auto x=values[i][a];if(!std::isfinite(x))throw std::runtime_error("Non-finite histogram weight");
      sums[i][a]+=x;for(unsigned b=0;b<channels;++b)products[i][a*channels+b]+=x*values[i][b];
    }
  }
};
struct Source {
  double mass,width,lo,hi,ecm,omega;
  std::pair<double,double> parent(double u)const {
    const double a=std::atan((lo*lo-mass*mass)/(mass*width));
    const double b=std::atan((hi*hi-mass*mass)/(mass*width));
    const double t=std::tan(a+(b-a)*u),s=mass*mass+mass*width*t;
    return {s,mass*width*(1+t*t)*(b-a)};
  }
};
std::vector<std::pair<double,double>> gauss(unsigned n){
  std::vector<std::pair<double,double>> rule;
  for(unsigned i=0;i<n;++i){double x=std::cos(M_PI*(i+.75)/(n+.5)),derivative=0;
    for(unsigned j=0;j<100;++j){double a=1,b=0;for(unsigned k=1;k<=n;++k){const auto old=a;a=((2*k-1)*x*a-(k-1)*b)/k;b=old;}
      derivative=n*(x*a-b)/(x*x-1);const auto dx=a/derivative;x-=dx;if(std::abs(dx)<1e-15)break;}
    rule.push_back({(1+x)/2,1/((1-x*x)*derivative*derivative)});
  }return rule;
}
// Independently angle-integrated SANC real amplitudes fix the inclusive
// normalization at EACH parent Q, without using any leptonff decay kernel.
struct InclusiveNormalization {
  int pid;double lo,step=1,omega;SANC::Process& proc;
  std::vector<std::array<double,2>> table;
  double quadrature_error=0,angular_error=0;
  double exact(double Q,bool up,unsigned nodes){
    double Z[8]{};Z[1]=Q*Q;Z[4]=up?.15:.05;
    double born=0,sv=0;proc.Born(Z,&born);proc.Virt(Z,&sv);
    if(!(born>0))throw std::runtime_error("Invalid inclusive SANC Born normalization");
    const double delta_sv=sv/born;
    Z[0]=.6;proc.Born(Z,&born);proc.Virt(Z,&sv);
    angular_error=std::max(angular_error,std::abs(sv/born-delta_sv));
    double integral_born=0,integral_sv=0;
    for(const auto& [u,w]:gauss(4)){Z[0]=2*u-1;proc.Born(Z,&born);proc.Virt(Z,&sv);
      integral_born+=2*w*born;integral_sv+=2*w*sv;}
    const double s=Q*Q,maximum=1-std::pow(proc.Mf3+proc.Mf4,2)/s,L=std::log(maximum/omega);
    long double real=0;
    for(const auto& [u,w]:gauss(nodes)){
      double shat=s,tau=omega*std::exp(L*u),spr=s*(1-tau),hard=0;
      if(pid==222)cc_ha_1413_1615_1spr_(&shat,&spr,&hard);
      else if(pid==212)cc_ha_1314_1516_1spr_(&shat,&spr,&hard);
      else if(up)nc_fsr_ha_1313_1616_1spr_(&shat,&spr,&hard);
      else nc_fsr_ha_1414_1616_1spr_(&shat,&spr,&hard);
      real+=w*s*tau*L*hard;
    }
    return (integral_sv+double(real))/integral_born;
  }
  InclusiveNormalization(int process,double lower,double upper,double cutoff,SANC::Process& p):
      pid(process),lo(lower),omega(cutoff),proc(p){
    const unsigned count=std::ceil((upper-lower)/step);step=(upper-lower)/count;
    for(unsigned i=0;i<=count;++i){const double Q=lo+i*step;
      const double up=exact(Q,true,512),down=pid==202?exact(Q,false,512):up;
      table.push_back({down,up});
      if(i==0 || i==count || i==count/2){quadrature_error=std::max(quadrature_error,std::abs(up-exact(Q,true,256)));}
    }
    if(quadrature_error>2e-6){
      std::cerr<<"Normalization angular variation "<<angular_error<<", quadrature shift "<<quadrature_error<<'\n';
      throw std::runtime_error("Inclusive SANC normalization quadrature needs refinement");}
  }
  double operator()(double Q,double flavour)const{
    const unsigned channel=int(flavour*10)+1;
    const bool up=(channel<6?channel%2==0:channel%2==1);
    const double x=(Q-lo)/step;const auto i=std::min(table.size()-2,std::size_t(x));const double f=x-i;
    return (1-f)*table[i][up]+f*table[i+1][up];
  }
};
}
int main(int argc,char** argv){try{
  if(argc<8)throw std::runtime_error("Usage: fsr_reference output.json pid trials seed omega mode(fo|showers|analytic|all) pdf [key=value ...]");
  if(std::ifstream(argv[1]))throw std::runtime_error("Refusing to overwrite output");
  const auto start=std::chrono::steady_clock::now();
  const int pid=std::stoi(argv[2]);const bool neutral=pid==202;
  if(pid!=202 && pid!=212 && pid!=222)throw std::runtime_error("Only muonic W-/W+/Z supported");
  const std::uint64_t trials=std::stoull(argv[3]),seed=std::stoull(argv[4]);
  if(trials<2)throw std::runtime_error("At least two trials required");
  const double omega=std::stod(argv[5]);const std::string mode=argv[6];
  const bool fo=mode=="fo"||mode=="all",showers=mode=="showers"||mode=="all",analytic=mode=="analytic"||mode=="all";
  if(!(fo||showers||analytic) || !(omega>0 && omega<.01))throw std::runtime_error("Invalid mode/soft cutoff");
  Inputs cfg;
  for(const auto& [k,v]:std::map<std::string,std::string>{{"pid",argv[2]},{"ecm","13000"},
      {"pdfname",argv[7]},{"iscale","1"},{"iborn","-1"},{"iqed","4"},{"iew","0"},
      {"iqcd","0"},{"ifgg","0"},{"iqedlli","0"},{"iqedllf","0"},{"iho","0"},
      {"gfscheme","0"},{"irun","0"},{"igmz","0"},{"mmo","0.1056583755"},
      {"alphai","137.035999084"},{"m34min",neutral?"50":"60"},{"m34max","200"}})cfg.v[k]=v;
  cfg.set("ome",omega);
  for(int i=8;i<argc;++i)cfg.set(std::string(argv[i]));
  if(cfg.get_l("iqed")!=4 || cfg.get_l("iew") || cfg.get_l("iqcd") || cfg.get_l("iborn")!=-1 ||
     cfg.get_l("iqedlli") || cfg.get_l("iqedllf") || cfg.get_l("iho") || cfg.get_l("gfscheme") ||
     cfg.get_l("iscale")!=1 || cfg.get_l("ifgg") || cfg.get_l("irun"))
    throw std::runtime_error("Unsupported physics flags: fixed-alpha pure FSR only");
  SANC::Info info(cfg);SANC::PartonDistributionsLHAPDF pdf;pdf.init(info.PDFname,info.PDFmember);
  SANC::Process proc(info,&pdf);SANC::ModulesInterface modules(info,proc);
  ffini_();
  Source source{neutral?info.MZ:info.MW,neutral?info.Wz:info.Ww,
                cfg.get_d("m34min"),cfg.get_d("m34max"),info.Ecm,omega};
  if(!(source.lo>2*info.Mmo && source.hi>source.lo && source.hi<source.ecm))throw std::runtime_error("Invalid parent interval");
  Histograms hist(neutral);
  // SANC orders W+ as (neutrino, antimuon), W- as (muon, antineutrino).
  const unsigned charged=pid==222?4:3,partner=pid==222?3:4;
  std::unique_ptr<InclusiveNormalization> normalization;
  if(fo)normalization.reset(new InclusiveNormalization(pid,source.lo,source.hi,omega,proc));
  const SpinDecay kind=neutral?SpinDecay::neutral_vector:pid==222?SpinDecay::wplus:SpinDecay::wminus;
  const double alpha=1/info.Alphai;
  if(showers && std::abs(info.Mmo-.1056583755)>1e-12)
    throw std::runtime_error("Shower comparison requires the common physical muon mass");
  std::unique_ptr<leptonff_dyturbo::PhotosDecay> photos;
  std::unique_ptr<leptonff_reference::PythiaDecay> pythia;
  const double shower_cutoff=cfg.v.count("shower_cutoff")?cfg.get_d("shower_cutoff"):1e-5;
  const bool mec=cfg.v.count("pythia_mec")?cfg.get_b("pythia_mec"):true;
  const bool quadrature=neutral && (mode=="analytic" || mode=="fo") && cfg.v.count("z_quadrature") && cfg.get_b("z_quadrature");
  if(showers){photos.reset(new leptonff_dyturbo::PhotosDecay(alpha,shower_cutoff,1e-7));
    pythia.reset(new leptonff_reference::PythiaDecay(alpha,shower_cutoff,mec));}
  const unsigned grid=cfg.v.count("ff_bins")?cfg.get_l("ff_bins"):65536;
  const double frag_factor=cfg.v.count("frag_scale_factor")?cfg.get_d("frag_scale_factor"):1;
  const double frag_scale=frag_factor*source.mass;
  if(!std::isfinite(frag_factor) || frag_factor<=0 || frag_scale<info.Mmo)
    throw std::runtime_error("Invalid fragmentation scale factor");
  const std::string mass_matching=cfg.v.count("mass_matching")?cfg.get_s("mass_matching"):"single_insertion";
  if(mass_matching!="single_insertion" && (mass_matching!="quadratic_finite" || !neutral || !analytic))
    throw std::runtime_error("Mass matching variation requires an analytic Z calculation");
  const auto matching=mass_matching=="quadratic_finite"?leptonff_dyturbo::MassMatching::quadratic_finite:
                                                     leptonff_dyturbo::MassMatching::single_insertion;
  std::unique_ptr<leptonff_dyturbo::MassRadiator> mass_ll,mass_nll;
  std::unique_ptr<DecayCalculator> ll,nll;
  std::unique_ptr<MeasurementPlan> plan;
  if(analytic){
    if(neutral){QEDParameters q;q.alpha=alpha;q.mu_frag=frag_scale;q.muon_mass=info.Mmo;q.mu_0=info.Mmo;
      LLGrid g;g.bins=grid;g.log_span=32;g.soft_log_cutoff=g.spacing()/2;
      mass_nll.reset(new leptonff_dyturbo::MassRadiator(q,g,PhotonicAccuracy::nll,matching));
      mass_ll.reset(new leptonff_dyturbo::MassRadiator(q,g,PhotonicAccuracy::ll,matching));
    }else{
      auto c=make_decay_configuration(kind,DecayPreset::photonic_ll,AlphaInput::on_shell(alpha),info.Mmo,frag_scale);
      c.grid.bins=grid;c.hard_samples=cfg.v.count("hard_samples")?cfg.get_l("hard_samples"):16;
      ll.reset(new DecayCalculator(c));c.mode=DecayMode::photonic_nll;nll.reset(new DecayCalculator(c));
      WLabCuts cuts;cuts.muon_pt_min=25;cuts.muon_eta_max=2.4;
      std::vector<Real> pt,eta,mt;for(unsigned i=0;i<=120;++i)pt.push_back(25+.25L*i);pt.push_back(INFINITY);
      for(unsigned i=0;i<=48;++i)eta.push_back(.05L*i);
      for(unsigned i=0;i<=200;++i)mt.push_back(.5L*i);mt.push_back(INFINITY);
      plan.reset(new MeasurementPlan(false,LabSpectrum::configured(false,{{"pt",pt},{"abseta",eta},{"mt",mt}},cuts)));
    }
  }
  std::mt19937_64 rng(seed);const auto uniform=[&](){return ((rng()>>12)+.5)/4503599627370496.;};
  double Z[8]{},H[8]{};std::array<SANC::Vec4,7> momenta;SANC::Vec4* p[7];for(unsigned i=0;i<7;++i)p[i]=&momenta[i];
  std::uint64_t negative_real=0,negative_hard=0,negative_pdf=0;double min_real=0;
  std::uint64_t quadrature_nodes=0;
  if(quadrature){
    const unsigned nq=cfg.v.count("mass_nodes")?cfg.get_l("mass_nodes"):16;
    const unsigned ny=cfg.v.count("rapidity_nodes")?cfg.get_l("rapidity_nodes"):32;
    if(nq<2 || nq>256 || ny<8 || ny>512)throw std::runtime_error("Invalid quadrature dimensions");
    const auto qr=gauss(nq),yr=gauss(ny),cr=gauss(3);
    const auto real_rule=gauss(cfg.v.count("real_nodes")?cfg.get_l("real_nodes"):8),tail_rule=gauss(128);
    const unsigned segments=std::ceil((source.hi-source.lo)/.125);
    const double h=(source.hi-source.lo)/segments;
    for(unsigned segment=0;segment<segments;++segment)for(const auto& [t,w]:qr){
      const double u=t*t*(3-2*t),M=source.lo+h*(segment+u),Y=std::log(info.Ecm/M);
      double density=0;std::array<double,2> luminosities{},amplitudes{};
      for(unsigned f=0;f<10;++f){
        double BornZ[8]{};BornZ[1]=M*M;BornZ[4]=(f+.5)/10;
        double angular=0;
        for(const auto& [v,cw]:cr){BornZ[0]=2*v-1;double B=0;proc.Born(BornZ,&B);angular+=2*cw*B;}
        double luminosity=0;
        for(const auto& [v,yw]:yr){BornZ[2]=(2*v-1)*Y;double lum=0;proc.PdfConvBORN(BornZ,&lum);luminosity+=2*Y*yw*lum;}
        density+=angular*luminosity;
        const bool up=(f<5?f%2==1:f%2==0);luminosities[up]+=luminosity;amplitudes[up]=angular;
      }
      const long double weight=2*M*density/(info.Ecm*info.Ecm)*h*6*t*(1-t)*w;
      const FourVector a{M/2,M/2,0,0},b{M/2,-M/2,0,0};const auto accepted=hist.measure(a,b);
      for(std::size_t i=0;i<hist.bins.size();++i)hist.sums[i][0]+=weight*accepted[i];
      if(analytic){
        hist.sums[0][5]+=weight;hist.sums[0][6]+=weight;hist.sums[0][7]+=weight;
        for(std::size_t i=2;i<hist.bins.size();++i){const auto& bin=hist.bins[i];
          const auto n=mass_nll->bin(M,bin.lo,bin.hi),l=mass_ll->bin(M,bin.lo,bin.hi);
          hist.sums[i][5]+=weight*n[1];hist.sums[i][6]+=weight*l[2];hist.sums[i][7]+=weight*n[2];
          hist.sums[1][5]+=weight*n[1];hist.sums[1][6]+=weight*l[2];hist.sums[1][7]+=weight*n[2];
        }
      }else{
        // For a mass-only measurement the unitary virtual contribution is
        // the endpoint subtraction of the independent SANC real spectrum.
        // Integrate the whole-plus action directly: no soft slicing bias.
        scales_.omega=0;
        for(std::size_t i=2;i<hist.bins.size();++i){const auto& bin=hist.bins[i];
          if(M<=bin.lo)continue;
          const bool born_bin=M<bin.hi;
          const double tau_lo=born_bin?1-bin.lo*bin.lo/(M*M):1-bin.hi*bin.hi/(M*M);
          const double tau_hi=born_bin?1-4*info.Mmo*info.Mmo/(M*M):1-bin.lo*bin.lo/(M*M);
          const auto& rule=born_bin?tail_rule:real_rule;
          const double log_ratio=std::log(tau_hi/tau_lo);
          long double delta=0;
          for(unsigned up=0;up<2;++up){long double integral=0;
            for(const auto& [v,rw]:rule){double s=M*M,tau=tau_lo*std::exp(v*log_ratio),spr=s*(1-tau),hard=0;
              if(up)nc_fsr_ha_1313_1616_1spr_(&s,&spr,&hard);
              else nc_fsr_ha_1414_1616_1spr_(&s,&spr,&hard);
              integral+=rw*s*tau*log_ratio*hard;}
            delta+=(born_bin?-1:1)*integral*luminosities[up];
          }
          const long double contribution=delta*(2*M/(info.Ecm*info.Ecm)*h*6*t*(1-t)*w);
          hist.sums[i][1]+=contribution;hist.sums[1][1]+=contribution;
        }
      }++quadrature_nodes;
    }
  }else for(std::uint64_t trial=0;trial<trials;++trial){
    const auto [s,jac_s]=source.parent(uniform());const double Q=std::sqrt(s),ymax=std::log(info.Ecm/Q),y=(2*uniform()-1)*ymax;
    const double c=2*uniform()-1,flavour=uniform(),azimuth=0;
    Z[0]=c;Z[1]=s;Z[2]=y;Z[3]=azimuth;Z[4]=flavour;
    double born=0,pdfs=0,sv=0;
    proc.PdfConvBORN(Z,&pdfs);proc.Born(Z,&born);
    const double jac_b=jac_s/(info.Ecm*info.Ecm)*(2*ymax)*2;
    const double weight=10*born*pdfs*jac_b;
    proc.GetFourMomentaBORN(Z,p);
    const auto a=read(*p[charged]),b=read(*p[partner]);const auto accepted=hist.measure(a,b);
    std::vector<Values> values(hist.bins.size());
    for(std::size_t i=0;i<values.size();++i)values[i][0]=weight*accepted[i];
    if(fo){
      proc.Virt(Z,&sv);const double wsv=10*sv*pdfs*jac_b;
      const double counter=weight*(*normalization)(Q,flavour);
      for(std::size_t i=0;i<values.size();++i){values[i][1]=(wsv-counter)*accepted[i];values[i][8]=counter*accepted[i];}
      const double tau_max=1-std::pow(proc.Mf3+proc.Mf4,2)/s;
      const double log_tau=std::log(tau_max/omega),tau=omega*std::exp(uniform()*log_tau),spr=s*(1-tau);
      const double beta=std::sqrt((spr-std::pow(proc.Mf3+proc.Mf4,2))*(spr-std::pow(proc.Mf3-proc.Mf4,2)))/(spr+info.Mmo*info.Mmo-(neutral?info.Mmo*info.Mmo:0));
      const double B=std::atanh(std::min(beta,1-1e-15)),v=B*(2*uniform()-1),cth4=std::tanh(v)/beta;
      const double jac_c=2*B/(beta*std::pow(std::cosh(v),2));
      H[0]=spr;H[1]=2*consts_.pi*uniform();H[2]=cth4;H[3]=2*uniform()-1;
      H[4]=s;H[5]=y;H[6]=0;H[7]=.25+.5*flavour;
      proc.PdfConvHARD(H,&pdfs);double hard=0;proc.Hard(H,&hard);
      const double jac_h=jac_s/(info.Ecm*info.Ecm)*(2*ymax)*s*tau*log_tau*(2*consts_.pi)*jac_c*2;
      const double wh=10*hard*pdfs*jac_h;
      if(hard<0)++negative_hard;if(pdfs<0)++negative_pdf;
      if(wh<0){++negative_real;min_real=std::min(min_real,wh);}
      proc.GetFourMomentaHARD(H,p);const auto f=hist.measure(read(*p[charged]),read(*p[partner]));
      for(std::size_t i=0;i<values.size();++i)values[i][2]=wh*f[i];
    }
    BornDecayPoint point{a,b,rng()};
    if(showers){
      const auto ph=photos->draw(kind,point,point.radiation_seed);
      const auto py=pythia->draw(kind,point,point.radiation_seed);
      const auto f=hist.measure(ph.first,ph.second),g=hist.measure(py.first,py.second);
      for(std::size_t i=0;i<values.size();++i){values[i][3]=weight*f[i];values[i][4]=weight*g[i];}
    }
    if(analytic){
      if(neutral){
        long double accounted1=0,accountedll=0,accountednll=0;
        for(std::size_t i=2;i<values.size();++i){const auto& bin=hist.bins[i];auto f=mass_nll->bin(Q,bin.lo,bin.hi);
          const auto l=mass_ll->bin(Q,bin.lo,bin.hi);
          values[i][5]=weight*f[1];values[i][6]=weight*l[2];values[i][7]=weight*f[2];
          accounted1+=f[1];accountedll+=l[2];accountednll+=f[2];}
        values[0][5]=values[0][6]=values[0][7]=weight;
        values[1][5]=weight*accounted1;values[1][6]=weight*accountedll;values[1][7]=weight*accountednll;
      }else{
        // The public operator expects massless Born parents. Keep the tiny
        // massive Born-map difference explicit through the shared baseline.
        const FourVector parent{a.energy+b.energy,a.px+b.px,a.py+b.py,a.pz+b.pz};
        const auto rest=rest_frame_boost(a,parent,-1);const Real norm=std::hypot(rest.pt(),rest.pz);
        auto first=rest_frame_boost({Q/2,Q*rest.px/(2*norm),Q*rest.py/(2*norm),Q*rest.pz/(2*norm)},parent,1);
        auto second=rest_frame_boost({Q/2,-Q*rest.px/(2*norm),-Q*rest.py/(2*norm),-Q*rest.pz/(2*norm)},parent,1);
        // Avoid cancellation in parent-first for a very backward daughter.
        first.energy=std::hypot(first.pt(),first.pz);second.energy=std::hypot(second.pt(),second.pz);
        const BornDecayPoint massless{first,second,point.radiation_seed};
        const auto pl=ll->evaluate(massless,*plan),pn=nll->evaluate(massless,*plan);
        values[0][5]=values[0][6]=values[0][7]=weight;
        values[1][5]=weight*pl.total.order1;values[1][6]=weight*pl.total.matched;values[1][7]=weight*pn.total.matched;
        for(std::size_t j=0;j<plan->spectrum().bins.size();++j){const auto& bin=plan->spectrum().bins[j];
          for(std::size_t i=2;i<values.size();++i)if(hist.bins[i].name==(bin.observable=="pt"?"ptl":bin.observable=="abseta"?"absetal":bin.observable) && std::abs(hist.bins[i].lo-bin.lower)<1e-8){
            values[i][5]=weight*pl.bins[j].order1;values[i][6]=weight*pl.bins[j].matched;values[i][7]=weight*pn.bins[j].matched;break;}}
      }
    }
    hist.add(values);
    if((trial+1)%100000==0)std::cerr<<"Completed "<<trial+1<<"/"<<trials<<'\n';
  }
  std::ofstream out(argv[1]);out<<std::setprecision(18);
  out<<"{\"schema\":\"leptonff-born-fsr-v1\",\"backend\":\"ReneSANCe-1.3.11/SANC\",\"mode\":\""<<mode
    <<"\",\"pid\":"<<pid<<",\"trials\":"<<trials<<",\"seed\":"<<seed<<",\"omega\":"<<omega
    <<",\"alpha\":"<<alpha<<",\"muon_mass\":"<<info.Mmo<<",\"parent_min\":"<<source.lo<<",\"parent_max\":"<<source.hi
    <<",\"pdf\":\""<<info.PDFname<<"\",\"pdf_member\":"<<info.PDFmember<<",\"ecm\":"<<info.Ecm
    <<",\"born_qcd\":true,\"weak\":false,\"isr\":false,\"ifi\":false,\"pairs\":false,\"running_alpha\":false"
    <<",\"shower_cutoff\":"<<shower_cutoff<<",\"pythia_mec\":"<<(mec?"true":"false")
    <<",\"photos_production_mec\":false,\"photos_legacy_w_correction\":true,\"photos_exponentiation_precision\":1e-7"
    <<",\"photos_version\":\""<<leptonff_dyturbo::PhotosDecay::version()<<"\",\"pythia_version\":\""<<leptonff_reference::PythiaDecay::version()
    <<"\",\"ff_bins\":"<<grid<<",\"frag_scale_factor\":"<<frag_factor<<",\"frag_scale\":"<<frag_scale
    <<",\"frag_scale_convention\":\"factor times pole mass; parent Q retained in hard matching\""
    <<",\"mass_matching\":\""<<mass_matching<<"\""
    <<",\"negative_real\":"<<negative_real<<",\"negative_hard\":"<<negative_hard<<",\"negative_pdf\":"<<negative_pdf<<",\"min_real\":"<<min_real
    <<",\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
    <<",\"integration\":\""<<(quadrature?"quadrature":"monte_carlo")<<"\",\"quadrature_nodes\":"<<quadrature_nodes
    <<",\"normalization\":\"unitary at each parent Q; independent SANC inclusive counterterm\""
    <<",\"normalization_quadrature_error\":"<<(normalization?normalization->quadrature_error:0)
    <<",\"normalization_angular_error\":"<<(normalization?normalization->angular_error:0)
    <<",\"settings\":{";
  bool first_setting=true;for(const auto& [k,v]:cfg.v){if(!first_setting)out<<',';first_setting=false;out<<std::quoted(k)<<':'<<std::quoted(v);}
  out<<"},\"channels\":[\"born\",\""<<(quadrature&&fo?"fo_plus_delta":"sv_delta")<<"\",\"real\",\"photos\",\"pythia\",\"our_order_alpha\",\"our_ll\",\"our_nll\",\"inclusive_counterterm\"],\"bins\":[";
  for(std::size_t i=0;i<hist.bins.size();++i){if(i)out<<',';const auto& bin=hist.bins[i];
    out<<"{\"observable\":\""<<bin.name<<"\",\"lower\":"<<bin.lo<<",\"upper\":"<<bin.hi<<",\"mean\":[";
    for(unsigned a=0;a<channels;++a){if(a)out<<',';out<<hist.sums[i][a]/(quadrature?1:trials);}
    out<<"],\"covariance\":[";for(unsigned a=0;a<channels;++a)for(unsigned b=0;b<channels;++b){if(a||b)out<<',';
      out<<(quadrature?0:(hist.products[i][a*channels+b]-hist.sums[i][a]*hist.sums[i][b]/trials)/(trials*(trials-1.L)));}
    out<<"]}";
  }out<<"]}\n";if(!out)throw std::runtime_error("Output write failed");
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}}
