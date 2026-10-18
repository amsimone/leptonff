// SPDX-License-Identifier: MIT
#pragma once
#include "leptonff/production_decay.hpp"
#include <string>
namespace leptonff_dyturbo {
struct PhotosDecay {
  static std::string version();
  PhotosDecay(double alpha,double cutoff,double exponentiation_precision=1e-4);
  std::pair<leptonff::FourVector,leptonff::FourVector> draw(
      leptonff::SpinDecay,const leptonff::BornDecayPoint&,std::uint64_t seed) const;
private:
  double cutoff_,precision_;
};
}
