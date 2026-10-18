#pragma once
#include "leptonff/spin_decay.hpp"
#include "leptonff/photonic_nll.hpp"

namespace leptonff {
// Fixed-alpha photonic reference with one finite insertion globally:
// U^n + a sum (U*d1) U^(n-1) + a H1[U^n]. Optional NLL adds
// a^2 L sum (U*P1_unitary) U^(n-1); at most one subleading insertion.
// The actual lattice O(alpha) overlap is replaced by continuum D1 and H1.
class PhotonicLabDressing {
public:
  PhotonicLabDressing(const QEDParameters&,const LLGrid&,Real minimum_z=1e-4L,
                     PhotonicAccuracy = PhotonicAccuracy::ll);
  const LLBoundaryNonSinglet& evolution() const { return evolution_; }
  const LabFragmentation& fragmentation() const { return fragmentation_; }
  const PhotonicNLL* nll_evolution() const { return nll_?&*nll_:nullptr; }
  // Pure LL measurement on a hard two- or three-body event, not U+a*d1.
  Real w_measurement(const FourVector&,const FourVector&,const WLabCuts& = {}) const;
  Real z_measurement(const FourVector&,const FourVector&,const FiducialCuts& = {}) const;
  // Coefficients multiplying a: H1[U^n]_lattice-H1[delta^n]_lattice.
  // Add these to the continuum bare projected coefficient, never to a
  // differently discretised overlap. Kernel scale uses the ORIGINAL hard M.
  Real w_projected_shift(const FourVector&,const FourVector&,const WLabCuts& = {}) const;
  Real z_projected_shift(const FourVector&,const FourVector&,const FiducialCuts& = {}) const;
  const std::vector<Real>& hard_source() const { return hard_; }
  const std::vector<Real>& w_extra_source() const { return w_extra_; }
  const std::vector<Real>& evolved_hard() const { return evolved_hard_; }
private:
  LLBoundaryNonSinglet evolution_;
  std::optional<PhotonicNLL> nll_;
  LabFragmentation fragmentation_;
  std::vector<Real> hard_,w_extra_,splitting_,evolved_hard_;
  SignedFractionTable u_,delta_,h_,hw_,p_,hu_,hwu_,pu_;
};
} // namespace leptonff
