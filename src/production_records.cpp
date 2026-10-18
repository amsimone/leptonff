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
GroupedIntegration::GroupedIntegration(std::size_t n):n_(n) {
  if(!n || n>8192)throw std::invalid_argument("Invalid replay channel count");
  sum_.resize(n);products_.resize(n*(n+1)/2);
}
void GroupedIntegration::add_group(const std::vector<Real>& x) {
  if(x.size()!=n_)throw std::invalid_argument("Replay channel mismatch");
  ++groups_;
  std::vector<std::size_t> active;
  for(std::size_t i=0;i<n_;++i) {if(!std::isfinite(x[i]))throw std::invalid_argument("Non-finite group response");sum_[i]+=x[i];if(x[i]!=0)active.push_back(i);}
  for(std::size_t j=0;j<active.size();++j)for(std::size_t i=0;i<=j;++i)products_[index(active[i],active[j])]+=x[active[i]]*x[active[j]];
}
void GroupedIntegration::end_block(const ProductionBlock& b) {
  if(b.trials<2 || b.stored_groups!=groups_ || groups_>b.trials || blocks_[{b.part,b.channel,b.iteration}])throw std::invalid_argument("Inconsistent/duplicate replay block");
  blocks_[{b.part,b.channel,b.iteration}]=true;
  auto& part=parts_[{b.part,b.channel}];
  if(!part.iterations) {part.sum.resize(n_);part.covariance.resize(products_.size());part.trials=b.trials;}
  // Equal iteration weights remain defined when the host resumes with a new
  // call count. Each block uses its OWN normalisation/zero-trial denominator.
  ++part.iterations;const Real calls=b.trials;
  for(std::size_t j=0;j<n_;++j) {
    part.sum[j]+=sum_[j];
    for(std::size_t i=0;i<=j;++i)part.covariance[index(i,j)]+=(calls*products_[index(i,j)]-sum_[i]*sum_[j])/(calls-1);
  }
  std::fill(sum_.begin(),sum_.end(),0);std::fill(products_.begin(),products_.end(),0);groups_=0;
}
IntegrationEstimate GroupedIntegration::estimate() const {
  if(groups_ || parts_.empty())throw std::logic_error("Incomplete/empty replay accumulation");
  IntegrationEstimate out;out.mean.resize(n_);out.covariance.resize(n_,std::vector<Real>(n_));
  for(const auto& entry:parts_) {
    const auto& p=entry.second;const Real n=p.iterations;
    for(std::size_t j=0;j<n_;++j) {
      out.mean[j]+=p.sum[j]/n;
      for(std::size_t i=0;i<=j;++i) {out.covariance[i][j]+=p.covariance[index(i,j)]/(n*n);if(i!=j)out.covariance[j][i]=out.covariance[i][j];}
    }
  }
  return out;
}
} // namespace leptonff
