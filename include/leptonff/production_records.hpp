#pragma once
#include "leptonff/production_decay.hpp"
#include <fstream>
#include <functional>
#include <map>
#include <string>
#include <tuple>

namespace leptonff {
struct ProductionMetadata {
  SpinDecay process=SpinDecay::wplus;
  std::string producer,physics_id,run_id,units="fb";
  std::uint64_t random_seed=0;
  bool full_born_angles=true,contains_qed=false;
  std::vector<std::string> weight_names{"central"};
  void validate() const;
};
struct ProductionTerm {BornDecayPoint point;std::vector<Real> weights;};
struct ProductionGroup {std::uint64_t id=0;std::vector<ProductionTerm> terms;};
struct ProductionBlock {
  int part=0,channel=0,iteration=0;
  std::uint64_t trials=0,stored_groups=0;
};

// Versioned, streaming, lossless text records. A group is written only after
// the host commits it; missing/zero trials remain in the BLOCK denominator.
// finish() is mandatory. An interrupted archive has no END and is rejected.
class ProductionWriter {
public:
  ProductionWriter(const std::string& filename,const ProductionMetadata&);
  void begin();
  void add(const ProductionTerm&);
  void commit(bool accepted=true);
  void block(const ProductionBlock&);
  void finish();
private:
  ProductionMetadata metadata_;
  std::ofstream file_;
  ProductionGroup pending_;
  std::uint64_t seen_=0,stored_=0;
  bool open_=false,finished_=false;
};
class ProductionReader {
public:
  explicit ProductionReader(const std::string& filename);
  const ProductionMetadata& metadata() const {return metadata_;}
  // Each callback receives a complete real/subtraction group. Returns the
  // footer for one block; false only after a verified END marker.
  bool next_block(ProductionBlock&,const std::function<void(const ProductionGroup&)>&);
private:
  ProductionMetadata metadata_;
  std::ifstream file_;
  bool ended_=false;
};

struct IntegrationEstimate {
  std::vector<Real> mean;
  std::vector<std::vector<Real>> covariance;
  std::vector<Real> variance;
  // Sparse/block layouts omit the dense matrix. Pairs include the diagonal;
  // absent entries are uncomputed, never implicitly zero covariance.
  std::vector<std::pair<std::size_t,std::size_t>> covariance_pairs;
  std::vector<Real> covariance_values;
};
enum class CovarianceMode { full,diagonal,blocks,selected };
struct CovarianceRequest {
  CovarianceMode mode=CovarianceMode::full;
  std::size_t block_size=0;
  std::vector<std::pair<std::size_t,std::size_t>> pairs;
};
// Sum counterevents before squaring. Host weights are already normalised
// integration contributions; trials are NOT another divisor of their sum.
class GroupedIntegration {
public:
  explicit GroupedIntegration(std::size_t channels,const CovarianceRequest& = {});
  void add_group(const std::vector<Real>&);
  void end_block(const ProductionBlock&);
  IntegrationEstimate estimate() const;
private:
  std::size_t n_;
  CovarianceRequest request_;
  std::vector<Real> sum_,products_;
  std::uint64_t groups_=0;
  struct Part {std::size_t iterations=0;std::uint64_t trials=0;std::vector<Real> sum,covariance;};
  std::map<std::pair<int,int>,Part> parts_;
  std::map<std::tuple<int,int,int>,bool> blocks_;
};
} // namespace leptonff
