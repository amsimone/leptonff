#include "leptonff/response_table.hpp"
#include "leptonff/decay_setup.hpp"
#include "parallel.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace leptonff {
namespace {
std::string identity(const DecayCalculator& c,const MeasurementPlan& m,const DecayAngularGrid& g) {
  std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(std::numeric_limits<Real>::max_digits10)
    <<decay_configuration_identity(c.configuration())<<" cartesian-lab-boost-v1 "<<g.cosine_nodes<<' '<<g.azimuth_nodes<<' '<<m.neutral()<<' ';
  for(const auto& v:m.spectrum().selections)s<<v.w.muon_pt_min<<' '<<v.w.muon_eta_max<<' '<<v.w.neutrino_pt_min<<' '<<v.w.transverse_mass_min<<' '
    <<v.z.pt_min<<' '<<v.z.eta_max<<' '<<v.z.mass_min<<' '<<v.z.mass_max<<';';
  s<<" bins ";for(const auto& b:m.spectrum().bins)s<<std::quoted(b.observable)<<' '<<b.lower<<' '<<b.upper<<' '<<b.positive<<' '<<b.negative<<';';
  return s.str();
}
void validate_node(const DecayResponseNode& n) {
  const auto& p=n.parent;
  for(Real v:{p.energy,p.px,p.py,p.pz})if(!std::isfinite(v))throw std::invalid_argument("Non-finite response node");
  if(p.energy<=0 || p.energy*p.energy<=p.px*p.px+p.py*p.py+p.pz*p.pz)throw std::invalid_argument("Response node must be timelike");
}
void dimensions(std::size_t nodes,std::size_t outputs) {
  if(!nodes || !outputs || outputs>100000 || nodes>10000000/(90*outputs))throw std::invalid_argument("Invalid/oversized response table");
}
// FNV-1a detects accidental corruption; it is not an authenticity signature.
void hash(std::uint64_t& h,const std::string& bytes) {for(unsigned char c:bytes){h^=c;h*=1099511628211ULL;}}
std::string record(const DecayResponseNode& node,const SpinDecayResponse& response) {
  std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(std::numeric_limits<Real>::max_digits10);
  const auto& p=node.parent;s<<"NODE "<<p.energy<<' '<<p.px<<' '<<p.py<<' '<<p.pz<<' '<<node.radiation_seed<<'\n';
  const auto put=[&](const auto& value){for(const auto& tensor:value)for(const auto& row:tensor)for(const auto& z:row)s<<z.real()<<' '<<z.imag()<<'\n';};
  put(response.total);for(const auto& b:response.bins)put(b);return s.str();
}
std::string header(const std::string& id,std::size_t n) {std::ostringstream s;s.imbue(std::locale::classic());s<<"LEPTONFF_RESPONSE 1\n"<<std::quoted(id)<<'\n'<<n<<'\n';return s.str();}
}
DecayResponseTable DecayResponseTable::build(const DecayCalculator& c,const MeasurementPlan& m,
    const std::vector<DecayResponseNode>& nodes,const DecayAngularGrid& grid,unsigned threads) {
  dimensions(nodes.size(),1+m.spectrum().bins.size());
  DecayResponseTable out;out.identity_=leptonff::identity(c,m,grid);out.nodes_=nodes;out.responses_.resize(nodes.size());
  for(const auto& n:nodes)validate_node(n);
  detail::parallel_for(nodes.size(),threads,[&](std::size_t i){out.responses_[i]=c.response(nodes[i].parent,m,nodes[i].radiation_seed,grid);});
  return out;
}
void DecayResponseTable::save(const std::string& filename) const {
  if(std::filesystem::exists(filename))throw std::invalid_argument("Response table output already exists");
  std::ofstream f(filename);f.imbue(std::locale::classic());if(!f)throw std::runtime_error("Cannot write response table");
  std::uint64_t h=14695981039346656037ULL;const auto head=header(identity_,nodes_.size());hash(h,head);f<<head;
  for(std::size_t i=0;i<nodes_.size();++i) {const auto data=record(nodes_[i],responses_[i]);hash(h,data);f<<data;}
  f<<"END "<<std::hex<<h<<'\n';f.close();if(!f)throw std::runtime_error("Response table write failed");
}
DecayResponseTable DecayResponseTable::load(const std::string& filename,const DecayCalculator& c,
    const MeasurementPlan& m,const DecayAngularGrid& grid) {
  std::ifstream f(filename);f.imbue(std::locale::classic());std::string magic,id;unsigned schema=0;std::size_t count=0;
  if(!(f>>magic>>schema) || magic!="LEPTONFF_RESPONSE" || schema!=1 || !(f>>std::quoted(id)>>count))throw std::invalid_argument("Invalid response table header");
  if(id!=leptonff::identity(c,m,grid))throw std::invalid_argument("Response table physics/measurement/numerical configuration mismatch");
  dimensions(count,1+m.spectrum().bins.size());
  DecayResponseTable out;out.identity_=id;out.nodes_.resize(count);out.responses_.resize(count);
  std::uint64_t h=14695981039346656037ULL;hash(h,header(id,count));
  for(std::size_t i=0;i<count;++i) {
    auto& n=out.nodes_[i];auto& r=out.responses_[i];r.bins.resize(m.spectrum().bins.size());
    if(!(f>>magic>>n.parent.energy>>n.parent.px>>n.parent.py>>n.parent.pz>>n.radiation_seed) || magic!="NODE")throw std::invalid_argument("Invalid response table node");
    validate_node(n);
    const auto get=[&](auto& value){
      for(auto& tensor:value) {
        for(auto& row:tensor)for(auto& z:row) {Real re=0,im=0;if(!(f>>re>>im) || !std::isfinite(re) || !std::isfinite(im))throw std::invalid_argument("Invalid/truncated response tensor");z={re,im};}
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)if(std::abs(tensor[a][b]-std::conj(tensor[b][a]))>1e-10L*std::max(Real{1},std::abs(tensor[a][b])))throw std::invalid_argument("Non-Hermitian response table tensor");
      }
    };
    get(r.total);for(auto& b:r.bins)get(b);hash(h,record(n,r));
  }
  std::uint64_t expected=0;
  if(!(f>>magic>>std::hex>>expected) || magic!="END" || expected!=h)throw std::invalid_argument("Incomplete/corrupt response table checksum");
  f>>std::ws;if(!f.eof())throw std::invalid_argument("Trailing response table data");return out;
}
} // namespace leptonff
