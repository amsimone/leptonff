#pragma once
#include "leptonff/mass_spectrum.hpp"

namespace leptonff {
enum class LeptonObservable { transverse_momentum, pseudorapidity };
// Shared numerical controls for a joint pT/eta calculation. Insert all
// one-lepton FF contacts into one panel partition, used by BOTH spectra.
CoupledDecayOptions lepton_spectrum_options(Real mass,const FiducialCuts& cuts,
    const std::vector<Real>& pt_edges,const std::vector<Real>& eta_edges,
    const CoupledDecayOptions& base = {});
// Both observed muons fill the histogram for EACH selected opposite-charge
// pair: sum = 2*pair rate, for every component and perturbative coefficient.
// Signed eta uses the exact reflection symmetry of THIS unpolarised-at-rest
// model. It must not be used unchanged for a boosted or polarised DY parent.
class CoupledLeptonSpectrum {
public:
  CoupledLeptonSpectrum(Real mass, const FiducialCuts& cuts, LeptonObservable observable,
      const std::vector<Real>& edges, const QEDXSpaceGrid& grid,
      const CoupledDecayOptions& options = {});
  // Reuse the bin/result layout of the mass spectrum. Here total is the PAIR
  // rate; sum contains two entries per pair. Closure is sum-2*total.
  CoupledMassSpectrumResult match(const QEDXSpaceEvolution& evolution,
      const IntegrationOptions& fixed_order_options = {}) const;
  const CoupledFiducialDecay& pair_operator() const { return cumulative_.front(); }
  const std::vector<CoupledFiducialDecay>& cumulative_operators() const { return cumulative_; }
private:
  LeptonObservable observable_;
  std::vector<Real> edges_, thresholds_;
  std::vector<CoupledFiducialDecay> cumulative_;
  // -1 means zero, -2 means the shared pair total; others index cumulative_.
  std::vector<int> indices_;
};
} // namespace leptonff
