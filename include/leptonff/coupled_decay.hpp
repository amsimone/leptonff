#pragma once
#include "leptonff/xspace_evolution.hpp"
#include "leptonff/fiducial_vector.hpp"

namespace leptonff {
struct CoupledDecayOptions {
  std::size_t ff_panels = 4; // four fixed cubic interpolation nodes per panel
  HardDecayOptions geometry;
  IntegrationOptions projection{1e-9L,1e-9L,8192};
  bool average_mass_cells = false; // optional conservative resolved-shell bin integration
  // Integrate every endpoint/resolved geometry entry on the original FF
  // lattice. Cubic compression is then used only for two resolved fractions.
  bool exact_endpoint_strips = false;
  // Additional panel boundaries in the observed FF fraction (0,1). Align
  // interpolation with measurement contacts while retaining a common linear
  // projection for all cumulative measurements.
  std::vector<Real> ff_fraction_breaks;
  void validate() const;
};

struct CoupledDecayMatch {
  Real born = 0, ff = 0, charge_flipped_pairs = 0;
  // Coefficients before multiplication by a_hard. Photon includes BOTH primary charges.
  Real muon_projected = 0, muon_geometry = 0;
  Real photon_projected = 0, photon_geometry = 0;
  Real a_initial = 0, a_hard = 0;
  Real lattice_coefficient1 = 0, continuum_coefficient1 = 0, coefficient2 = 0;
  Real overlap_correction = 0, raw = 0, matched = 0, fixed_order = 0;
  Real compression_moment_residual = 0, compression_absolute_norm = 0;
};

// Reusable hard operator for symmetric cuts on an unpolarised vector at rest.
// All O(alpha) hard-parent channels; inclusive distinct-parent opposite-charge
// pairs, no extra-pair veto. Same explicit MSbar input for FFs and hard a(M).
// Positive mass_min is required. Finite geometry uses a LINEAR signed cubic
// projection of truncated unweighted FFs; endpoint atoms remain explicit.
// Hard mass projections/convolutions are causal on physical support, no FFT.
class CoupledFiducialDecay {
public:
  CoupledFiducialDecay(Real mass, const FiducialCuts& cuts,
      const QEDXSpaceGrid& grid, const CoupledDecayOptions& options = {});
  CoupledDecayMatch match(const QEDXSpaceEvolution& evolution,
      const IntegrationOptions& fixed_order_options = {}) const;
  // Restrict the mass window while RETAINING this operator's FF support,
  // hard source projections and linear interpolation nodes. This is useful
  // for additive cumulative histograms. Other cuts and numerical controls stay
  // fixed; the requested window must lie inside the original window.
  CoupledFiducialDecay restricted_mass_window(Real lower, Real upper) const;
  // Mean one-lepton cumulative, keeping the original pair cuts and the SAME
  // FF projection/compression. This replaces any prior additional selection.
  CoupledFiducialDecay with_lepton_selection(const SingleLeptonSelection& selection) const;
  Real minimum_fraction() const { return minimum_; }
  const SingleLeptonSelection& lepton_selection() const { return selection_; }
  const std::vector<Real>& fractions() const { return fractions_; }
  const std::vector<Real>& ff_panel_edges() const { return panel_edges_; }
  // Prefix arrays represent z*H(z) dz, at mu_frag=M. Their endpoint actions
  // include the entire plus subtraction, including the unrepresented tail.
  const std::vector<Real>& muon_projection() const { return hm_; }
  const std::vector<Real>& photon_projection() const { return hg_; }
  const std::vector<std::vector<Real>>& muon_geometry() const { return gm_; }
  const std::vector<std::vector<Real>>& photon_geometry() const { return gg_; }
  // Direct G(z_j,1)/z_j and (for photons) G(1,z_j)/z_j, including j=0.
  // Populated when exact_endpoint_strips is enabled; source arrays are z*D.
  const std::vector<Real>& muon_endpoint_geometry() const { return em_; }
  const std::vector<Real>& photon_first_endpoint_geometry() const { return eg_first_; }
  const std::vector<Real>& photon_second_endpoint_geometry() const { return eg_second_; }
  // Exposed for reproducible numerical audits (same layout as fractions).
  std::vector<Real> compress(const std::vector<Real>& weighted) const;
private:
  void build_geometry();
  Real mass_, minimum_;
  FiducialCuts cuts_;
  QEDXSpaceGrid grid_;
  CoupledDecayOptions options_;
  SingleLeptonSelection selection_;
  std::size_t prefix_;
  std::vector<Real> fractions_, panel_edges_, hm_, hg_, pm_, pg_;
  std::vector<Real> em_, eg_first_, eg_second_;
  std::vector<std::vector<Real>> compression_, gm_, gg_;
};
} // namespace leptonff
