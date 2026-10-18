#pragma once
#include "leptonff/lab_fragmentation.hpp"
#include <string>

namespace leptonff {
// A physical histogram is a linear combination of nested fiducial actions.
// All terms must use the SAME hard/radiation draws before any squaring.
struct LabSelection {
  WLabCuts w;
  FiducialCuts z;
  Real accepts(bool neutral,const FourVector&,const FourVector&) const;
};
struct LabSpectrumBin {
  std::string observable;
  Real lower,upper; // +infinity denotes the retained overflow bin
  int positive,negative; // selection indices; -1 denotes the zero functional
};
struct ObservableBinning {
  std::string observable;
  std::vector<Real> edges;
};
struct LabSpectrum {
  std::vector<LabSelection> selections; // index zero is the original rate
  std::vector<LabSpectrumBin> bins;
  static LabSpectrum reference(bool neutral,const WLabCuts& = {},const FiducialCuts& = {});
  // Each requested observable partitions its entire fiducial support. Infinite
  // overflow is retained for pT/mT; duplicate/unsupported observables rejected.
  static LabSpectrum configured(bool neutral,const std::vector<ObservableBinning>&,
                               const WLabCuts& = {},const FiducialCuts& = {});
  static std::vector<ObservableBinning> reference_binning(bool neutral,unsigned subdivision=1);
  static std::vector<ObservableBinning> read_binning(const std::string& filename);
  Real bin(std::size_t,const std::vector<Real>& cumulative) const;
};

// Immutable validated measurement, reusable with different physics contexts.
// Selections differing only in pseudorapidity share the same fraction action.
class MeasurementPlan {
public:
  MeasurementPlan(bool neutral,LabSpectrum);
  bool neutral() const {return neutral_;}
  const LabSpectrum& spectrum() const {return spectrum_;}
  const std::vector<std::size_t>& representatives() const {return representatives_;}
  bool accepts_eta(std::size_t,const FourVector&,const FourVector&) const;
private:
  bool neutral_;
  LabSpectrum spectrum_;
  std::vector<std::size_t> representatives_;
};
} // namespace leptonff
