#include "leptonff/lab_fragmentation.hpp"
#include "leptonff/spin_decay.hpp"
#include "leptonff/photonic_dressing.hpp"
#include "leptonff/coupled_lab.hpp"
#include "leptonff/lab_spectrum.hpp"
#include "leptonff/hadronic_input.hpp"
#include "leptonff/production_decay.hpp"
#include "leptonff/production_records.hpp"
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <atomic>
#include <random>
#include <cstring>
#include <iomanip>
#include <string>
#ifdef LEPTONFF_WITH_PHOTOS
#include "photos_comparison.hpp"
#endif

namespace {
using namespace leptonff;
Real setting(const char* name,Real fallback) {
  const char* value=std::getenv(name);
  if(!value)return fallback;
  std::size_t used=0;const Real parsed=std::stold(value,&used);
  if(used!=std::string(value).size() || !std::isfinite(parsed))
    throw std::invalid_argument(std::string("Invalid setting ")+name);
  return parsed;
}
struct Engine {
  WLabCuts wcuts;
  FiducialCuts zcuts;
  std::unique_ptr<DecayCalculator> calculator;
  const LabFragmentation* ff=nullptr;
  const PhotonicLabDressing* dressing=nullptr;
  const CoupledLabDressing* coupled=nullptr;
  Real mu_frag=0,soft=0,collinear=0;
  int hard_samples=0;
  bool spectra=false;
  explicit Engine(int process) {
    wcuts.muon_pt_min=setting("LEPTONFF_PT_MIN",25);
    wcuts.muon_eta_max=setting("LEPTONFF_ETA_MAX",2.4);
    wcuts.neutrino_pt_min=setting("LEPTONFF_NU_PT_MIN",0);
    wcuts.transverse_mass_min=setting("LEPTONFF_MT_MIN",0);
    zcuts.pt_min=wcuts.muon_pt_min;zcuts.eta_max=wcuts.muon_eta_max;
    zcuts.mass_min=setting("LEPTONFF_MASS_MIN",60);
    zcuts.mass_max=setting("LEPTONFF_MASS_MAX",120);
    wcuts.validate();zcuts.validate();
#ifdef LEPTONFF_WITH_PHOTOS
    if(wcuts.muon_pt_min!=25 || wcuts.muon_eta_max!=Real(2.4) ||
       wcuts.neutrino_pt_min!=0 || wcuts.transverse_mass_min!=0 ||
       zcuts.mass_min!=60 || zcuts.mass_max!=120)
      throw std::invalid_argument("PHOTOS pilot uses fixed pT25/eta2.4, no neutrino/mT cuts, Z mass60-120");
#endif
    MomentEvolutionParameters p;
    p.mu_frag=setting("LEPTONFF_MU_FRAG",process==31?91.1876L:80.385L);
    p.alpha_msbar_at_mu0=setting("LEPTONFF_ALPHA",p.alpha_msbar_at_mu0);
    const char* flavour_value=std::getenv("LEPTONFF_FLAVOURS");
    const std::string flavour=flavour_value?flavour_value:"one-family";
    const Real threshold_factor=setting("LEPTONFF_THRESHOLD_FACTOR",1);
    HadronicAlphaInput hadronic;
    if(flavour=="hadronic-pointlike") {
      if(threshold_factor!=1)throw std::invalid_argument("Use hadronic matching scale for this profile");
      const char* file=std::getenv("LEPTONFF_HADRONIC_FILE");if(file)hadronic=HadronicAlphaInput::from_file(file);
      hadronic.alpha_s_reference=setting("LEPTONFF_ALPHA_S_REFERENCE",hadronic.alpha_s_reference);
      hadronic.delta_had5_os+=setting("LEPTONFF_HADRONIC_SHIFT",0)*hadronic.uncertainty;
      configure_hadronic_pointlike(p,hadronic,p.alpha_msbar_at_mu0,setting("LEPTONFF_HADRONIC_SCALE",2));
      if(setting("LEPTONFF_QCD_OFF",0)==1)p.qcd->alpha_s_reference=0;
    }else if(flavour=="leptons")p.flavour_scheme=QEDFlavourScheme::leptons(threshold_factor);
    else if(flavour=="custom") {
      const char* file=std::getenv("LEPTONFF_FLAVOUR_FILE");
      if(!file)throw std::invalid_argument("Custom flavours need a flavour file");
      if(threshold_factor!=1)throw std::invalid_argument("Custom threshold scales must be specified in the file");
      p.flavour_scheme=QEDFlavourScheme::from_file(file);
    }else if(flavour!="one-family" || threshold_factor!=1)throw std::invalid_argument("Invalid flavour profile/threshold factor");
    const char* scheme_value=std::getenv("LEPTONFF_ALPHA_SCHEME");
    const std::string input_scheme=scheme_value?scheme_value:"msbar";
    const Real input_alpha=setting("LEPTONFF_ALPHA",1/137.035999084L);
    if(flavour=="hadronic-pointlike") {
      if(input_scheme!="on-shell")throw std::invalid_argument("HVP anchor requires on-shell input");
    }else if(input_scheme=="on-shell") {
      QEDFlavourScheme reference;reference.fermions={{"mu",-1,p.muon_mass,1,0}};
      p.alpha_msbar_at_mu0=leptonic_alpha_msbar(input_alpha,p.mu_0,p.flavour_scheme?*p.flavour_scheme:reference);
    }else if(input_scheme!="msbar")throw std::invalid_argument("Invalid alpha input scheme");
    mu_frag=p.mu_frag;
    const Real samples=setting("LEPTONFF_HARD_SAMPLES",0);
    if(samples!=std::floor(samples) || samples<0 || samples>4096)
      throw std::invalid_argument("Invalid hard sample count");
    hard_samples=static_cast<int>(samples);
    soft=setting("LEPTONFF_HARD_SOFT",1e-5L);
    collinear=setting("LEPTONFF_HARD_COLLINEAR",1e-5L);
    if(soft<=0 || soft>=.01L || collinear<=0 || collinear>=.01L)
      throw std::invalid_argument("Invalid finite-geometry cutoffs");
    QEDXSpaceGrid grid;grid.log_span=setting("LEPTONFF_LOG_SPAN",p.qcd?64:40);
    const Real bins=setting("LEPTONFF_GRID_BINS",4096);
    if(bins!=std::floor(bins) || bins<16 || bins>262144)throw std::invalid_argument("Invalid FF grid bins");
    grid.bins=static_cast<std::size_t>(bins);
    const bool photonic=setting("LEPTONFF_PHOTONIC_LL",0)==1;
    const bool photonic_nll=setting("LEPTONFF_PHOTONIC_NLL",0)==1;
    if(photonic_nll && (!photonic || setting("LEPTONFF_DRESS_HARD",0)!=1 || hard_samples==0))
      throw std::invalid_argument("Photonic NLL requires photonic evolution and finite hard dressing");
    const bool coupled_hard=setting("LEPTONFF_COUPLED_HARD",0)==1;
    if(p.flavour_scheme && !coupled_hard)throw std::invalid_argument("Variable flavours require coupled hard matching in this adapter");
    if(photonic && input_scheme!="msbar")throw std::invalid_argument("Photonic mode uses its explicit fixed alpha input, no conversion");
    if(coupled_hard && (photonic || hard_samples==0))
      throw std::invalid_argument("Coupled hard matching needs coupled evolution and hard samples");
    spectra=setting("LEPTONFF_SPECTRA",0)==1;
    if(spectra && !coupled_hard && !(photonic && setting("LEPTONFF_DRESS_HARD",0)==1 && hard_samples>0))
      throw std::invalid_argument("Spectra require coupled hard or photonic dressed hard mode");
    const bool fraction_cells=setting("LEPTONFF_FRACTION_CELLS",0)==1;
    if(fraction_cells && !coupled_hard)throw std::invalid_argument("Conservative cells require coupled hard mode");
    DecayConfiguration config;
    config.process=process==31?SpinDecay::neutral_vector:process==1?SpinDecay::wplus:SpinDecay::wminus;
    config.mode=photonic?(photonic_nll?DecayMode::photonic_nll:DecayMode::photonic_ll):
      coupled_hard?DecayMode::coupled_nll:DecayMode::ff_only;
    config.evolution=p;config.grid=grid;config.hard_samples=hard_samples;
    config.soft_cutoff=soft;config.collinear_cutoff=collinear;
    config.dress_photonic_hard=setting("LEPTONFF_DRESS_HARD",0)==1;
    config.fraction_measure=fraction_cells?FractionMeasure::conservative_cells:FractionMeasure::atoms;
    calculator=std::make_unique<DecayCalculator>(config);
    ff=&calculator->fragmentation();dressing=calculator->photonic();coupled=calculator->coupled();
    if(!photonic && setting("LEPTONFF_DRESS_HARD",0)!=0)
      throw std::invalid_argument("Hard dressing pilot requires the fixed-alpha photonic branch");
    std::cerr<<"leptonff MCFM: Born-parent FF action; grid="
      <<grid.bins<<", log_span="<<grid.log_span<<", mu_frag="<<p.mu_frag<<", alpha(mu0)="<<p.alpha_msbar_at_mu0
      <<"; mode="<<(photonic?"fixed-alpha photonic LL + continuum D1":p.flavour_scheme?"variable-flavour MSbar":"coupled one-family MSbar")<<'\n';
    if(hard_samples)std::cerr<<"leptonff finite decay: unitary O(alpha), spin tomography, additive H1; samples="
      <<hard_samples<<", soft="<<soft<<", collinear="<<collinear
      <<"; Z charge-symmetric rates only; hard dressing absent\n";
    if(coupled)std::cerr<<"leptonff coupled hard: "<<(p.flavour_scheme?"variable-flavour pure-QED":"one-family")<<" MSbar NLL, all hard parents, inclusive pairs; common alpha(M), continuum overlap; not full SM/hadronic content or an aligned PHOTOS comparison\n";
    if(p.qcd)std::cerr<<std::setprecision(18)<<"leptonff hadronic pointlike: matching_scale="<<p.hadronic_matching->scale
      <<"; extra_coefficient="<<p.hadronic_matching->coefficient<<"; delta_had5_os="<<hadronic.delta_had5_os
      <<"; input_uncertainty="<<hadronic.uncertainty<<"; alpha_s_reference="<<p.qcd->alpha_s_reference
      <<"; reference_scale="<<hadronic.reference_scale<<"; explicit pointlike boundary, LO QCD, no mixed kernels or fitted low-scale FF\n";
    if(p.flavour_scheme) {
      const auto run=run_flavour_coupling(p,p.mu_frag);
      std::cerr<<std::setprecision(18)<<"leptonff flavour input: profile="<<flavour
        <<"; alpha_input="<<input_alpha<<"; scheme="<<input_scheme<<"; alpha_msbar_mu0="<<p.alpha_msbar_at_mu0
        <<"; alpha_msbar_mu_frag="<<2*std::acos(-1.L)*run.a_final
        <<"; beta0_integral_mu_frag="<<run.beta0_integral<<"; finite time-like threshold matching enabled\n";
      for(const auto& f:p.flavour_scheme->fermions)std::cerr<<"leptonff fermion: "<<f.name<<", charge="<<f.charge
        <<", colours="<<f.colours<<", mass="<<f.mass<<", threshold="<<f.threshold()<<", pointlike="<<f.pointlike<<'\n';
    }
    if(fraction_cells)std::cerr<<"leptonff fraction cells: signed resolved uniform cells; number/energy conserved, exact endpoint atom, exact pair cut integrals; continuum overlap and C2 use the same measure\n";
    if(spectra)std::cerr<<"leptonff spectra: physical observed W pT/mT/absolute eta; Z pair mass/minimum muon pT/maximum absolute eta; shared currents and radiation draws\n";
    if(dressing)std::cerr<<"leptonff photonic dressing: U^n + one dressed boundary or H1 globally; continuum O(alpha) restored; "
      <<(photonic_nll?"unitary photonic P1 insertion; no pairs/running":"no P1/pairs")<<'\n';
    if(photonic_nll) {
      const auto& d=dressing->nll_evolution()->diagnostics();
      std::cerr<<std::setprecision(18)<<"leptonff photonic NLL: convention=unitary-conditioned; source_N1="
        <<d.source_number_moment<<"; raw_inclusive_factor="<<d.raw_inclusive_factor
        <<"; kernel_number_residual="<<d.kernel_number_residual<<"; evolved_number_residual="<<d.evolved_number_residual
        <<"; wrap_bound="<<d.wrap_bound<<"; tail_bound="<<d.omitted_tail_bound<<'\n';
    }
  }
};
const Engine& engine_for(int process) {
  if(process!=1 && process!=6 && process!=31)throw std::invalid_argument("Unsupported MCFM process");
  static const int initialized_process=process;
  if(process!=initialized_process)throw std::invalid_argument("One process per adapter executable is required");
  static const Engine engine(process);return engine;
}
struct ShellDiagnostics {
  std::atomic<double> relative{0},absolute{0};
  static void update(std::atomic<double>& target,double x) {
    double old=target.load();
    while(x>old && !target.compare_exchange_weak(old,x)){}
  }
  ~ShellDiagnostics() {
    std::cerr<<"leptonff MCFM shell projection: max_relative_energy_shift="<<relative.load()
      <<", max_absolute_energy_shift_GeV="<<absolute.load()<<'\n';
  }
} shell;
FourVector vector(const double* p) {
  const Real norm=std::hypot(std::hypot(Real(p[0]),Real(p[1])),Real(p[2]));
  const Real delta=std::abs(norm-p[3]);
  if(!std::isfinite(norm) || !std::isfinite(p[3]) || p[3]<=0 || delta>1e-6L*p[3])
    throw std::invalid_argument("MCFM parent exceeds massless mapping-roundoff tolerance");
  // MCFM dipole mappings can leave tiny off-shell residuals on soft leptons.
  // Preserve all three-momenta; repair only energy, with explicit diagnostics.
  ShellDiagnostics::update(shell.relative,double(delta/p[3]));
  ShellDiagnostics::update(shell.absolute,double(delta));
  return {norm,p[0],p[1],p[2]};
}
std::uint64_t hard_seed(const double* a,const double* b) {
  std::uint64_t value=0x7a2f9c4d1103b6e5ULL;
  for(const auto* p:{a,b})for(int k=0;k<4;++k) {
    std::uint64_t bits;std::memcpy(&bits,p+k,sizeof(bits));
    value^=bits;value*=1099511628211ULL;
  }
  return value;
}
std::array<Real,4> hard_geometry(const Engine& e,int process,const FourVector& a,const FourVector& b,
                   const double* ref3,const double* ref4) {
  if(!e.hard_samples || e.ff->a_initial()==0)return {};
  const FourVector q{a.energy+b.energy,a.px+b.px,a.py+b.py,a.pz+b.pz};
  const Real M=invariant_mass(a,b),pi=std::acos(-1.L);
  const auto rest=rest_frame_boost(a,q,-1);
  const Real norm=std::hypot(rest.pt(),rest.pz);
  const auto kind=process==31?SpinDecay::neutral_vector:process==1?SpinDecay::wplus:SpinDecay::wminus;
  const auto density=spin_angular_projector(kind,{rest.px/norm,rest.py/norm,rest.pz/norm});
  const auto measurement=[&](const FourVector& mu,const FourVector& other)->Real {
    if(process==31)return e.zcuts.accepts(mu,other)?1:0;
    return mu.pt()>e.wcuts.muon_pt_min && std::abs(mu.eta())<e.wcuts.muon_eta_max &&
      (e.wcuts.neutrino_pt_min==0 || other.pt()>e.wcuts.neutrino_pt_min) &&
      (e.wcuts.transverse_mass_min==0 || transverse_mass(mu,other)>e.wcuts.transverse_mass_min)?1:0;
  };
  const auto evolved_measurement=[&](const FourVector& mu,const FourVector& other)->Real {
    return process==31?e.dressing->z_measurement(mu,other,e.zcuts)
                      :e.dressing->w_measurement(mu,other,e.wcuts);
  };
  std::mt19937_64 rng(hard_seed(ref3,ref4));
  const auto uniform=[&](){return Real((rng()>>11)+.5L)/9007199254740992.L;};
  const Real edge=std::acos(1-e.collinear),theta_span=pi-2*edge;
  const Real tmin=std::sqrt(e.soft),tspan=std::sqrt(1-1e-8L)-tmin;
  std::array<Real,4> sum{};
  for(int k=0;k<e.hard_samples;++k) {
    const Real t=tmin+tspan*uniform(),z=1-t*t;
    const Real theta=edge+theta_span*uniform(),c=std::cos(theta);
    // Keep draw order explicit across compilers and mapped counterevents.
    const Real orientation=2*uniform()-1,phi=2*pi*uniform(),twist=2*pi*uniform();
    const auto rotation=decay_rotation(orientation,phi,twist);
    const Real pair_mass=M*std::sqrt(z);
    const Real reference=process!=31 || (pair_mass>e.zcuts.mass_min && pair_mass<e.zcuts.mass_max)?1:0;
    if(e.coupled) {
      const auto measurements=[&](const FourVector& mu,const FourVector& other,const FourVector& photon) {
        return std::array<Real,3>{measurement(mu,other),
          process==31?e.coupled->z_real(mu,other,photon,e.zcuts):e.coupled->w_real(mu,other,photon,e.wcuts),
          process==31?e.coupled->z_real(mu,other,photon,e.zcuts,true):e.coupled->w_real(mu,other,photon,e.wcuts,true)};
      };
      const auto g=decay_hard_geometry_channels(kind,M,z,c,rotation,q,measurements,reference);
      const Real jacobian=2*t*tspan*theta_span*std::sin(theta);
      sum[0]+=jacobian*spin_contract(density,g[0]);
      sum[2]+=jacobian*spin_contract(density,g[1]);
      sum[3]+=jacobian*spin_contract(density,g[2]);
      continue;
    }
    const auto g=e.dressing?decay_hard_geometry_pair(kind,M,z,c,rotation,q,measurement,evolved_measurement,reference)
        :std::array<SpinTensor,2>{decay_hard_geometry(kind,M,z,c,rotation,q,measurement,reference),{}};
    const Real jacobian=2*t*tspan*theta_span*std::sin(theta);
    for(int i=0;i<2;++i)sum[i]+=jacobian*spin_contract(density,g[i]);
  }
  for(auto& value:sum)value/=e.hard_samples;
  return sum;
}
}

namespace {
struct ArchiveState {
  std::unique_ptr<ProductionWriter> writer;
  int process=0;
  ArchiveState() {
    const char* filename=std::getenv("LEPTONFF_PRODUCTION_ARCHIVE");
    if(!filename || !*filename)return;
    const char* threads=std::getenv("OMP_NUM_THREADS");
    if(!threads || std::string(threads)!="1")throw std::invalid_argument("Production archive requires one thread");
    const auto required=[](const char* key){const char* value=std::getenv(key);if(!value || !*value)throw std::invalid_argument(std::string("Missing archive metadata: ")+key);return std::string(value);};
    process=std::stoi(required("LEPTONFF_ARCHIVE_PROCESS"));
    if(process!=1 && process!=6 && process!=31)throw std::invalid_argument("Invalid archive process");
    ProductionMetadata meta;meta.process=process==31?SpinDecay::neutral_vector:process==1?SpinDecay::wplus:SpinDecay::wminus;
    meta.producer="MCFM-10.3.1";meta.physics_id=required("LEPTONFF_ARCHIVE_PHYSICS_ID");meta.run_id=required("LEPTONFF_ARCHIVE_RUN_ID");
    meta.random_seed=std::stoull(required("LEPTONFF_ARCHIVE_RANDOM_SEED"));
    if(required("LEPTONFF_ARCHIVE_FULL_ANGLES")!="1")throw std::invalid_argument("Archive requires full Born angular coverage");
    writer=std::make_unique<ProductionWriter>(filename,meta);
    std::cerr<<"leptonff production archive: grouped pre-FSR points, central PDF/scale weight only; units fb\n";
  }
  ~ArchiveState() {if(writer)try{writer->finish();}catch(const std::exception& e){std::cerr<<"Incomplete production archive: "<<e.what()<<'\n';}}
};
ArchiveState& archive() {static ArchiveState state;return state;}
template<class F> void archive_call(F&& action) {
  try{auto& state=archive();if(state.writer)action(state);}
  catch(const std::exception& e){std::cerr<<"leptonff production archive failed: "<<e.what()<<'\n';std::abort();}
}
}
extern "C" void leptonff_archive_begin_() {archive_call([](ArchiveState& s){s.writer->begin();});}
extern "C" void leptonff_archive_end_(const int* accepted) {archive_call([&](ArchiveState& s){s.writer->commit(*accepted!=0);});}
extern "C" void leptonff_archive_add_(const int* process,const double* p3,const double* p4,const double* ref3,const double* ref4,const double* weight) {
  archive_call([&](ArchiveState& s){
    if(*process!=s.process)throw std::invalid_argument("Changed archive process");
    const BornDecayPoint point{map_massless_parent(*process==1?p4:p3),map_massless_parent(*process==1?p3:p4),radiation_seed(ref3,ref4)};
    s.writer->add({point,{*weight}});
  });
}
extern "C" void leptonff_archive_flush_(const int* stage,const int* part,const int* ips,const int* iteration,const std::int64_t* calls) {
  archive_call([&](ArchiveState& s){if(*stage<1)throw std::invalid_argument("Warm-up archive export is unsupported");s.writer->block({*part,*ips,*iteration,std::uint64_t(*calls),0});});
}

extern "C" int leptonff_mcfm_evaluate(int process,const double* p3,const double* p4,
    const double* ref3,const double* ref4,double* out) noexcept {
  try {
    if(process!=1 && process!=6 && process!=31)
      throw std::invalid_argument("Adapter supports only MCFM W+/W-/Z processes 1/6/31");
    // Each executable run has exactly one selected process. C++ static
    // initialization is synchronized; subsequent event evaluation is read-only.
    const auto& engine=engine_for(process);
    const auto a=vector(p3),b=vector(p4);
    const auto result=process==31?engine.ff->z(a,b,engine.zcuts)
      :engine.ff->w(process==1?b:a,process==1?a:b,engine.wcuts);
    out[0]=result.born;
    out[1]=result.fixed_order;
    out[2]=result.resummed;
    out[3]=result.matched;
    out[4]=result.matched-result.born;
    out[5]=result.matched-result.fixed_order;
    out[6]=engine.ff->a_initial()*engine.ff->a_initial()*result.coefficient2;
    for(int k=7;k<33;++k)out[k]=0;
    if(engine.hard_samples) {
      const auto charged=process==1?b:a,other=process==1?a:b;
      out[12]=engine.ff->a_initial()*(process==31?z_projected_decay_hard(a,b,engine.mu_frag,engine.zcuts)
                    :w_projected_decay_hard(charged,other,engine.mu_frag,engine.wcuts));
      const auto geometry=hard_geometry(engine,process,charged,other,ref3,ref4);
      out[13]=engine.ff->a_initial()*geometry[0];
      out[14]=out[12]+out[13];out[15]=result.fixed_order+out[14];out[16]=result.matched+out[14];
      if(engine.coupled) {
        const auto h=process==31?engine.coupled->z_projected(a,b,engine.zcuts)
          :engine.coupled->w_projected(charged,other,engine.wcuts);
        const auto matched=CoupledLabDressing::match(result,h,{geometry[0],geometry[2],geometry[3]});
        out[26]=h.a_initial*h.continuum+h.a_hard*h.evolved-h.a_initial*h.lattice;
        out[27]=h.a_hard*geometry[2];out[28]=out[26]+out[27];
        out[29]=result.matched+out[28];out[30]=h.a_initial*h.a_initial*matched.coefficient2;
        out[31]=out[29]-out[15];
      }
      if(engine.dressing) {
        const auto boundary=process==31?engine.dressing->fragmentation().z(a,b,engine.zcuts)
                    :engine.dressing->fragmentation().w(charged,other,engine.wcuts);
        out[19]=boundary.matched-result.matched;
        if(engine.dressing->nll_evolution())out[6]=engine.ff->a_initial()*engine.ff->a_initial()*boundary.coefficient2;
        out[20]=engine.ff->a_initial()*(process==31?engine.dressing->z_projected_shift(a,b,engine.zcuts)
                    :engine.dressing->w_projected_shift(charged,other,engine.wcuts));
        out[21]=engine.ff->a_initial()*geometry[1];
        out[22]=out[19]+out[20]+out[21];out[23]=out[16]+out[22];out[25]=out[23]-out[15];
      }
    }
#ifdef LEPTONFF_WITH_PHOTOS
    double photos[4];
    leptonff_photos_comparison(process,p3,p4,ref3,ref4,photos);
    out[7]=photos[0];out[8]=photos[1];out[9]=photos[2];
    out[10]=photos[1]-result.matched;out[11]=photos[3];
    if(engine.hard_samples)out[17]=photos[1]-out[16];
    if(engine.hard_samples)out[18]=photos[1]-out[15];
    if(engine.coupled)out[32]=photos[1]-out[29];
    if(engine.dressing)out[24]=photos[1]-out[23];
#else
    (void)ref3;(void)ref4;
#endif
    return 0;
  }catch(const std::exception& e) {
    std::cerr<<"leptonff MCFM adapter failed: "<<e.what()<<'\n';return 1;
  }
}

// Eight columns per physical bin, Fortran layout (8,17). Rate ABI is unchanged.
extern "C" int leptonff_mcfm_evaluate_all(int process,const double* p3,const double* p4,
    const double* ref3,const double* ref4,double* rates,double* spectra) noexcept {
  if(leptonff_mcfm_evaluate(process,p3,p4,ref3,ref4,rates))return 1;
  try {
    std::fill(spectra,spectra+8*17,0);
    const auto& e=engine_for(process);if(!e.spectra)return 0;
    const bool neutral=process==31;
    const auto definition=LabSpectrum::reference(neutral,e.wcuts,e.zcuts);
    const auto a=vector(process==1?p4:p3),b=vector(process==1?p3:p4);
    const auto kind=neutral?SpinDecay::neutral_vector:process==1?SpinDecay::wplus:SpinDecay::wminus;
    const auto predictions=e.calculator->cumulative({a,b,hard_seed(ref3,ref4)},definition.selections);
    const auto n=definition.selections.size();
    std::vector<std::array<Real,8>> values(n);
#ifdef LEPTONFF_WITH_PHOTOS
    std::vector<DecayMeasurement> bare;
    for(const auto& selection:definition.selections)
      bare.push_back([=](const FourVector& u,const FourVector& v){return selection.accepts(neutral,u,v);});
    const auto photos=leptonff_photos_batch(process,p3,p4,ref3,ref4,bare);
#endif
    for(std::size_t i=0;i<n;++i) {
      const auto& p=predictions[i];
      values[i]={p.born,p.order1,p.matched,p.order2,p.higher,0,0,0};
#ifdef LEPTONFF_WITH_PHOTOS
      values[i][5]=photos[i][0];values[i][6]=photos[i][1];values[i][7]=photos[i][1]-p.matched;
#endif
    }
    // The same samples and selections must reproduce the original rate ABI.
    for(int j=0;j<8;++j) {
      const std::array<int,8> indices=e.dressing?std::array<int,8>{0,15,23,6,25,7,8,24}:std::array<int,8>{0,15,29,30,31,7,8,32};
      if(std::abs(values[0][j]-rates[indices[j]])>1e-9L*(1+std::abs(values[0][j])))
        throw std::runtime_error("Differential/base-rate identity failed");
    }
    if(definition.bins.size()!=17)throw std::runtime_error("Unexpected spectrum size");
    for(std::size_t k=0;k<17;++k)for(int j=0;j<8;++j) {
      const auto& b=definition.bins[k];
      spectra[8*k+j]=double((b.positive<0?0:values[b.positive][j])-(b.negative<0?0:values[b.negative][j]));
    }
    return 0;
  }catch(const std::exception& ex) {
    std::cerr<<"leptonff spectra failed: "<<ex.what()<<'\n';return 1;
  }
}
