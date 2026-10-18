#pragma once
#include "leptonff/spin_decay.hpp"

namespace leptonff {
// Same cluster projection as projected_decay_hard_kernel, but photon fraction
// x=1-z. REGULAR (no plus/delta); small-x number integral does not exist.
Real projected_decay_photon_kernel(SpinDecay,Real x,Real M,Real mu_frag);

struct CoupledLabProjection {
  Real a_initial=0,a_hard=0,beta0_log_hard=0;
  Real continuum=0,lattice=0,evolved=0,first=0;
};
struct CoupledLabMatch {
  Real raw=0,overlap=0,matched=0,fixed_order=0,coefficient2=0;
};
// Geometry entries: bare, full coupled measurement, coefficient of a_initial
// of that measurement. All use ONE shared real current and collinear map.
using CoupledGeometry = std::array<SpinTensor,3>;
using ParentMeasurements = std::function<std::array<Real,3>(
    const FourVector&,const FourVector&,const FourVector&)>;
CoupledGeometry decay_hard_geometry_channels(SpinDecay,Real M,Real z,Real c,
    const SpatialRotation&,const FourVector& parent,const ParentMeasurements&,
    Real orientation_independent_reference=0);

// Multiple measurements share the current, boost, rotation and sector maps.
std::vector<CoupledGeometry> decay_hard_geometry_channels_batch(SpinDecay,Real M,Real z,Real c,
    const SpatialRotation&,const FourVector& parent,const std::vector<ParentMeasurements>&,
    const std::vector<Real>& orientation_independent_references);

// General linear channels with a separate orientation-independent control
// for each channel (e.g. bare and evolved-minus-bare in photonic matching).
std::vector<CoupledGeometry> decay_hard_geometry_linear_batch(SpinDecay,Real M,Real z,Real c,
    const SpatialRotation&,const FourVector& parent,const std::vector<ParentMeasurements>&,
    const std::vector<std::array<Real,3>>& orientation_independent_references);

// One-family MSbar NLL reference: full P0/P1, running, d1 and every hard-parent
// channel. Inclusive tagged W leptons / distinct-parent opposite-charge Z
// pairs. Positive pT and Z pair-mass cuts exclude unresolved same-parent pairs.
// Retained higher orders include products of finite insertions; no NNLO claim.
class CoupledLabDressing {
public:
  CoupledLabDressing(SpinDecay,const MomentEvolutionParameters& = {},
                    const QEDXSpaceGrid& = {},Real minimum_z=1e-4L,
                    FractionMeasure = FractionMeasure::atoms);
  const QEDXSpaceEvolution& evolution() const { return evolution_; }
  const LabFragmentation& fragmentation() const { return fragmentation_; }
  Real a_hard(Real original_mass) const;
  CoupledLabProjection w_projected(const FourVector&,const FourVector&,const WLabCuts& = {}) const;
  CoupledLabProjection z_projected(const FourVector&,const FourVector&,const FiducialCuts& = {}) const;
  // Three-body hard measurement, including photon conversion and wrong-sign
  // primaries. first_order selects its explicit O(a_initial) coefficient.
  Real w_real(const FourVector& mu,const FourVector& nu,const FourVector& photon,
              const WLabCuts& = {},bool first_order=false) const;
  Real z_real(const FourVector& minus,const FourVector& plus,const FourVector& photon,
              const FiducialCuts& = {},bool first_order=false) const;
  static CoupledLabMatch match(const LabFFResult&,const CoupledLabProjection&,
                              const std::array<Real,3>& contracted_geometry);
private:
  enum Table { identity,same,opposite,gamma,c1,cg,h0,p0,hs,ho,hg,ps,po,pg,hc,hcg,pc,pcg,charged,count };
  SpinDecay kind_;
  QEDXSpaceEvolution evolution_;
  LabFragmentation fragmentation_;
  std::vector<SignedFractionTable> t_;
  CoupledLabProjection projection(Real M,const std::function<Real(Table,Table)>& action,bool pair) const;
};
} // namespace leptonff
