// SPDX-License-Identifier: MIT
#pragma once
#include "leptonff/production_decay.hpp"
#include <memory>
namespace leptonff_reference {
class PythiaDecay {
public:
  PythiaDecay(double alpha, double cutoff, bool matrix_element);
  ~PythiaDecay();
  std::pair<leptonff::FourVector,leptonff::FourVector> draw(
      leptonff::SpinDecay,const leptonff::BornDecayPoint&,std::uint64_t);
  static std::string version();
private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};
}
