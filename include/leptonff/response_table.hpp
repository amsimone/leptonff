#pragma once
#include "leptonff/production_decay.hpp"
#include <string>

namespace leptonff {
struct DecayResponseNode {FourVector parent;std::uint64_t radiation_seed=0;};
// Exact-node cache, independent of production weights. It introduces no
// interpolation approximation. Angular and hard-sampling convergence remain
// independent requirements of the underlying response calculation.
class DecayResponseTable {
public:
  static DecayResponseTable build(const DecayCalculator&,const MeasurementPlan&,
      const std::vector<DecayResponseNode>&,const DecayAngularGrid& = {},unsigned threads=1);
  static DecayResponseTable load(const std::string& filename,const DecayCalculator&,
      const MeasurementPlan&,const DecayAngularGrid& = {});
  void save(const std::string& filename) const;
  std::size_t size() const {return nodes_.size();}
  const DecayResponseNode& node(std::size_t index) const {return nodes_.at(index);}
  const SpinDecayResponse& response(std::size_t index) const {return responses_.at(index);}
  const std::string& identity() const {return identity_;}
private:
  DecayResponseTable()=default;
  std::string identity_;
  std::vector<DecayResponseNode> nodes_;
  std::vector<SpinDecayResponse> responses_;
};
} // namespace leptonff
