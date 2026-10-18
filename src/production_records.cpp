#include "leptonff/production_records.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <set>
#include <stdexcept>

namespace leptonff {
namespace {
int code(SpinDecay p) {return p==SpinDecay::wplus?1:p==SpinDecay::wminus?6:31;}
SpinDecay process(int p) {if(p!=1 && p!=6 && p!=31)throw std::invalid_argument("Invalid archive process");return p==1?SpinDecay::wplus:p==6?SpinDecay::wminus:SpinDecay::neutral_vector;}
void term(const ProductionTerm& p,std::size_t weights) {
  if(p.weights.size()!=weights)throw std::invalid_argument("Archive weight schema mismatch");
  for(Real w:p.weights)if(!std::isfinite(w))throw std::invalid_argument("Non-finite production weight");
  for(const auto& v:{p.point.first,p.point.second}) {
    for(Real x:{v.energy,v.px,v.py,v.pz})if(!std::isfinite(x))throw std::invalid_argument("Non-finite production momentum");
    if(v.energy<=0 || std::abs(v.energy*v.energy-v.px*v.px-v.py*v.py-v.pz*v.pz)>1e-9L*v.energy*v.energy)
      throw std::invalid_argument("Archive requires massless pre-FSR Born leptons");
  }
}
void write_vector(std::ostream& out,const FourVector& p) {out<<' '<<p.energy<<' '<<p.px<<' '<<p.py<<' '<<p.pz;}
void read_vector(std::istream& in,FourVector& p) {in>>p.energy>>p.px>>p.py>>p.pz;}
std::size_t index(std::size_t i,std::size_t j) {if(i>j)std::swap(i,j);return j*(j+1)/2+i;}
}
void ProductionMetadata::validate() const {
  if(process!=SpinDecay::wplus && process!=SpinDecay::wminus && process!=SpinDecay::neutral_vector)throw std::invalid_argument("Invalid archive decay");
  if(producer.empty() || physics_id.empty() || run_id.empty() || (units!="fb" && units!="pb"))throw std::invalid_argument("Incomplete production metadata");
  if(!full_born_angles || contains_qed)throw std::invalid_argument("Replay requires full Born angles and production without QED FSR");
  std::set<std::string> names;
  if(weight_names.empty() || weight_names.size()>1024)throw std::invalid_argument("Invalid production weight count");
  for(const auto& name:weight_names)if(name.empty() || !names.insert(name).second)throw std::invalid_argument("Invalid production weight names");
}
ProductionWriter::ProductionWriter(const std::string& name,const ProductionMetadata& m):metadata_(m) {
  m.validate();if(std::filesystem::exists(name))throw std::invalid_argument("Refusing to overwrite production archive: "+name);
  file_.open(name);if(!file_)throw std::runtime_error("Cannot create production archive");
  file_<<std::setprecision(std::numeric_limits<Real>::max_digits10);
  file_<<"LEPTONFF_PRODUCTION 1\nMETA "<<code(m.process)<<' '<<std::quoted(m.producer)<<' '<<std::quoted(m.physics_id)<<' '
    <<std::quoted(m.run_id)<<' '<<std::quoted(m.units)<<' '<<m.random_seed<<' '<<m.full_born_angles<<' '<<m.contains_qed<<' '<<m.weight_names.size()<<"\nWEIGHTS";
  for(const auto& name:m.weight_names)file_<<' '<<std::quoted(name);
  file_<<'\n';
}
void ProductionWriter::begin() {
  if(open_ || finished_)throw std::logic_error("Invalid archive group begin");
  open_=true;pending_={++seen_,{}};
}
void ProductionWriter::add(const ProductionTerm& p) {
  if(!open_)throw std::logic_error("Archive term outside integration group");
  term(p,metadata_.weight_names.size());
  if(pending_.terms.size()>=4096)throw std::runtime_error("Too many terms in subtraction group");
  if(!pending_.terms.empty() && p.point.radiation_seed!=pending_.terms.front().point.radiation_seed)
    throw std::invalid_argument("Subtraction terms must share the original radiation seed");
  pending_.terms.push_back(p);
}
void ProductionWriter::commit(bool accepted) {
  if(!open_)throw std::logic_error("Archive commit without begin");
  open_=false;
  if(accepted && !pending_.terms.empty()) {
    file_<<"GROUP "<<pending_.id<<' '<<pending_.terms.size()<<'\n';
    for(const auto& p:pending_.terms) {
      file_<<"POINT "<<p.point.radiation_seed;write_vector(file_,p.point.first);write_vector(file_,p.point.second);
      for(Real w:p.weights)file_<<' '<<w;
      file_<<'\n';
    }
    ++stored_;
  }
  pending_.terms.clear();
}
void ProductionWriter::block(const ProductionBlock& b) {
  if(open_ || finished_ || b.trials<2 || seen_>b.trials || b.iteration<1)throw std::logic_error("Invalid archive block footer");
  file_<<"BLOCK "<<b.part<<' '<<b.channel<<' '<<b.iteration<<' '<<b.trials<<' '<<stored_<<'\n';
  file_.flush();if(!file_)throw std::runtime_error("Production archive write failed");seen_=stored_=0;
}
void ProductionWriter::finish() {
  if(finished_)return;
  if(open_ || seen_ || stored_)throw std::logic_error("Uncommitted production archive data");
  file_<<"END\n";file_.flush();if(!file_)throw std::runtime_error("Production archive close failed");finished_=true;
}
ProductionReader::ProductionReader(const std::string& name):file_(name) {
  std::string token;int version=0,p=0;std::size_t weights=0;
  if(!(file_>>token>>version) || token!="LEPTONFF_PRODUCTION" || version!=1)throw std::invalid_argument("Unsupported production archive");
  if(!(file_>>token>>p>>std::quoted(metadata_.producer)>>std::quoted(metadata_.physics_id)>>std::quoted(metadata_.run_id)
       >>std::quoted(metadata_.units)>>metadata_.random_seed>>metadata_.full_born_angles>>metadata_.contains_qed>>weights)
      || token!="META" || !weights || weights>1024)throw std::invalid_argument("Malformed production metadata");
  metadata_.process=process(p);metadata_.weight_names.resize(weights);
  if(!(file_>>token) || token!="WEIGHTS")throw std::invalid_argument("Missing weight schema");
  for(auto& name:metadata_.weight_names)file_>>std::quoted(name);
  if(!file_)throw std::invalid_argument("Truncated production metadata");
  metadata_.validate();
}
bool ProductionReader::next_block(ProductionBlock& block,const std::function<void(const ProductionGroup&)>& consume) {
  if(ended_)return false;
  std::uint64_t groups=0,previous=0;std::string token;
  while(file_>>token) {
    if(token=="END") {
      if(groups)throw std::invalid_argument("Archive ends inside a block");
      if(file_>>token)throw std::invalid_argument("Trailing archive data");
      ended_=true;return false;
    }
    if(token=="BLOCK") {
      if(!(file_>>block.part>>block.channel>>block.iteration>>block.trials>>block.stored_groups) || block.trials<2 ||
          block.stored_groups!=groups || previous>block.trials || block.iteration<1)throw std::invalid_argument("Inconsistent archive block");
      return true;
    }
    if(token!="GROUP")throw std::invalid_argument("Unknown production archive record");
    ProductionGroup group;std::size_t count=0;
    if(!(file_>>group.id>>count) || group.id<=previous || !count || count>4096)throw std::invalid_argument("Invalid archive group");
    previous=group.id;group.terms.resize(count);
    for(auto& p:group.terms) {
      if(!(file_>>token>>p.point.radiation_seed) || token!="POINT")throw std::invalid_argument("Missing archive point");
      read_vector(file_,p.point.first);read_vector(file_,p.point.second);p.weights.resize(metadata_.weight_names.size());
      for(auto& w:p.weights)file_>>w;
      if(!file_)throw std::invalid_argument("Truncated production point");
      term(p,p.weights.size());
      if(p.point.radiation_seed!=group.terms.front().point.radiation_seed)throw std::invalid_argument("Inconsistent group radiation seeds");
    }
    consume(group);++groups;
  }
  throw std::invalid_argument("Interrupted production archive: missing END");
}
namespace {
template<class F> void covariance_pairs(std::size_t n,const CovarianceRequest& r,const F& f) {
  std::size_t k=0;
  if(r.mode==CovarianceMode::selected) {for(auto p:r.pairs)f(p.first,p.second,k++);return;}
  const auto size=r.mode==CovarianceMode::full?n:r.mode==CovarianceMode::diagonal?1:r.block_size;
  for(std::size_t start=0;start<n;start+=size)
    for(std::size_t j=start;j<std::min(n,start+size);++j)
      for(std::size_t i=start;i<=j;++i)f(i,j,k++);
}
}
GroupedIntegration::GroupedIntegration(std::size_t n,const CovarianceRequest& request):n_(n),request_(request) {
  if(!n || n>1000000)throw std::invalid_argument("Invalid replay channel count");
  if(request_.mode!=CovarianceMode::full && request_.mode!=CovarianceMode::diagonal &&
     request_.mode!=CovarianceMode::blocks && request_.mode!=CovarianceMode::selected)throw std::invalid_argument("Invalid covariance mode");
  if(request_.mode==CovarianceMode::blocks && (!request_.block_size || request_.block_size>n))throw std::invalid_argument("Invalid covariance block size");
  if(request_.mode==CovarianceMode::full && n>8192)throw std::invalid_argument("Full covariance exceeds 8192 channels; choose blocks, selected pairs or diagonal");
  if(request_.mode==CovarianceMode::selected) {
    for(auto& p:request_.pairs) {if(p.first>=n || p.second>=n)throw std::invalid_argument("Invalid selected covariance channel");if(p.first>p.second)std::swap(p.first,p.second);}
    for(std::size_t i=0;i<n;++i)request_.pairs.emplace_back(i,i);
    std::sort(request_.pairs.begin(),request_.pairs.end());
    request_.pairs.erase(std::unique(request_.pairs.begin(),request_.pairs.end()),request_.pairs.end());
  }
  std::size_t count=0;
  if(request_.mode==CovarianceMode::selected)count=request_.pairs.size();
  else {
    const auto b=request_.mode==CovarianceMode::full?n:request_.mode==CovarianceMode::diagonal?1:request_.block_size;
    count=(n/b)*(b*(b+1)/2)+(n%b)*((n%b)+1)/2;
  }
  if(count>100000000)throw std::invalid_argument("Requested covariance exceeds 100 million stored entries");
  sum_.resize(n);products_.resize(count);
}
void GroupedIntegration::add_group(const std::vector<Real>& x) {
  if(x.size()!=n_)throw std::invalid_argument("Replay channel mismatch");
  for(auto v:x)if(!std::isfinite(v))throw std::invalid_argument("Non-finite group response");
  ++groups_;
  for(std::size_t i=0;i<n_;++i)sum_[i]+=x[i];
  if(request_.mode==CovarianceMode::selected) {
    for(std::size_t k=0;k<request_.pairs.size();++k) {auto p=request_.pairs[k];products_[k]+=x[p.first]*x[p.second];}
    return;
  }
  const auto b=request_.mode==CovarianceMode::full?n_:request_.mode==CovarianceMode::diagonal?1:request_.block_size;
  std::vector<std::size_t> active;active.reserve(b);
  for(std::size_t start=0;start<n_;start+=b) {
    active.clear();
    for(std::size_t i=start;i<std::min(n_,start+b);++i)if(x[i]!=0)active.push_back(i);
    const auto offset=(start/b)*(b*(b+1)/2);
    for(std::size_t j=0;j<active.size();++j)for(std::size_t i=0;i<=j;++i)
      products_[offset+index(active[i]-start,active[j]-start)]+=x[active[i]]*x[active[j]];
  }
}
void GroupedIntegration::end_block(const ProductionBlock& b) {
  if(b.trials<2 || b.stored_groups!=groups_ || groups_>b.trials || blocks_[{b.part,b.channel,b.iteration}])throw std::invalid_argument("Inconsistent/duplicate replay block");
  blocks_[{b.part,b.channel,b.iteration}]=true;
  auto& part=parts_[{b.part,b.channel}];
  if(!part.iterations) {part.sum.resize(n_);part.covariance.resize(products_.size());part.trials=b.trials;}
  ++part.iterations;const Real calls=b.trials;
  for(std::size_t i=0;i<n_;++i)part.sum[i]+=sum_[i];
  covariance_pairs(n_,request_,[&](std::size_t i,std::size_t j,std::size_t k){part.covariance[k]+=(calls*products_[k]-sum_[i]*sum_[j])/(calls-1);});
  std::fill(sum_.begin(),sum_.end(),0);std::fill(products_.begin(),products_.end(),0);groups_=0;
}
IntegrationEstimate GroupedIntegration::estimate() const {
  if(groups_ || parts_.empty())throw std::logic_error("Incomplete/empty replay accumulation");
  IntegrationEstimate out;out.mean.resize(n_);out.variance.resize(n_);
  const bool full=request_.mode==CovarianceMode::full;
  if(full)out.covariance.resize(n_,std::vector<Real>(n_));
  else {
    out.covariance_pairs.reserve(products_.size());out.covariance_values.resize(products_.size());
    covariance_pairs(n_,request_,[&](std::size_t i,std::size_t j,std::size_t){out.covariance_pairs.emplace_back(i,j);});
  }
  for(const auto& entry:parts_) {
    const auto& p=entry.second;const Real n=p.iterations;
    for(std::size_t j=0;j<n_;++j)out.mean[j]+=p.sum[j]/n;
    covariance_pairs(n_,request_,[&](std::size_t i,std::size_t j,std::size_t k){
      const Real v=p.covariance[k]/(n*n);
      if(i==j)out.variance[i]+=v;
      if(full) {out.covariance[i][j]+=v;if(i!=j)out.covariance[j][i]=out.covariance[i][j];}
      else out.covariance_values[k]+=v;
    });
  }
  return out;
}
} // namespace leptonff
