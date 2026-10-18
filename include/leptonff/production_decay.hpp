#pragma once
#include "leptonff/coupled_lab.hpp"
#include "leptonff/lab_spectrum.hpp"
#include "leptonff/photonic_dressing.hpp"
#include <cstdint>
#include <memory>

namespace leptonff {
enum class DecayMode { ff_only, photonic_ll, photonic_nll, coupled_nll };
struct DecayConfiguration {
  SpinDecay process=SpinDecay::wplus;
  DecayMode mode=DecayMode::photonic_nll;
  MomentEvolutionParameters evolution;
  QEDXSpaceGrid grid;
  unsigned hard_samples=16;
  Real soft_cutoff=1e-5L,collinear_cutoff=1e-5L,minimum_z=1e-4L;
  bool dress_photonic_hard=true;
  FractionMeasure fraction_measure=FractionMeasure::atoms;
  void validate() const;
};

// first = charged W lepton, or negative Z lepton; second = neutrino or
// positive Z lepton. Host weights live outside this decay functional.
// A subtraction group MUST use a common seed across its mapped points.
struct BornDecayPoint {
  FourVector first,second;
  std::uint64_t radiation_seed=0;
};
std::uint64_t radiation_seed(const double* reference3,const double* reference4);
FourVector map_massless_parent(const double* px_py_pz_energy,Real relative_tolerance=1e-6L);

struct DecayPrediction {
  Real born=0,order1=0,matched=0,order2=0,higher=0;
  std::array<Real,5> values() const {return {born,order1,matched,order2,higher};}
};
struct DecaySpectrumPrediction {
  DecayPrediction total;
  std::vector<DecayPrediction> bins;
};

// rho and response use rest-frame Cartesian axes defined by a rotation-free
// boost from the LAB axes. W: Hermitian 3x3; Z: parity-even real symmetric
// tensors only. A Collins-Soper/helicity host must rotate its rho first.
struct SpinDecayResponse {
  std::array<SpinTensor,5> total;
  std::vector<std::array<SpinTensor,5>> bins;
  DecaySpectrumPrediction contract(const SpinTensor& production) const;
};
struct DecayAngularGrid {unsigned cosine_nodes=16,azimuth_nodes=32;};

// Immutable, generator-independent, and safe for concurrent evaluation.
// No environment variables, host globals, PHOTOS, or production PDF dependence.
class DecayCalculator {
public:
  explicit DecayCalculator(const DecayConfiguration&);
  const DecayConfiguration& configuration() const;
  const LabFragmentation& fragmentation() const;
  const PhotonicLabDressing* photonic() const;
  const CoupledLabDressing* coupled() const;
  std::vector<DecayPrediction> cumulative(const BornDecayPoint&,
                                         const std::vector<LabSelection>&) const;
  DecaySpectrumPrediction evaluate(const BornDecayPoint&,const LabSpectrum&) const;
  // Direct production-tensor route. Angular quadrature errors and finite
  // geometry sampling are separate controls; refine both for precision.
  SpinDecayResponse response(const FourVector& parent,const LabSpectrum&,
                            std::uint64_t seed,const DecayAngularGrid& = {}) const;
private:
  struct Impl;
  std::shared_ptr<const Impl> impl_;
};
} // namespace leptonff
