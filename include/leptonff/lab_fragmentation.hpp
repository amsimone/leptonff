#pragma once
#include "leptonff/xspace_evolution.hpp"
#include "leptonff/decay.hpp"
#include "leptonff/evolution.hpp"
#include <optional>

namespace leptonff {
class PhotonicNLL;
class FractionPairWorkspace;
struct WLabCuts {
  Real muon_pt_min = 25, muon_eta_max = 2.4L;
  Real neutrino_pt_min = 0, transverse_mass_min = 0;
  void validate() const;
};
// Conventional pT-based mT, also usable with massive bare PHOTOS muons.
Real transverse_mass(const FourVector& muon,const FourVector& neutrino);

struct LabFFResult {
  Real born = 0, resummed = 0, lattice_coefficient1 = 0;
  Real continuum_coefficient1 = 0, coefficient2 = 0;
  Real fixed_order = 0, overlap_correction = 0, matched = 0;
};

enum class FractionMeasure { atoms, conservative_cells };

// Signed number-weight prefixes on the original lattice, including its atom.
// Bounds are STRICT. Default nodes are atoms. Conservative cells spread each
// resolved weight uniformly around its node, preserving number and energy;
// the z=1 atom is always exact. No clipping or renormalisation. Queries
// below minimum_z are rejected: the full coupled number density is unbounded.
class SignedFractionTable {
public:
  SignedFractionTable(Real spacing,const std::vector<Real>& number_weights,
                      Real minimum_z,FractionMeasure = FractionMeasure::atoms);
  Real interval(Real lower,Real upper = 2) const;
  // Integral with z1>lower1, z2>lower2, product_min<z1*z2<product_max.
  // Exact for the chosen atomic/uniform-cell measures, including both atoms.
  Real pair(const SignedFractionTable& second,Real lower1,Real lower2,
            Real product_min,Real product_max) const;
  Real inverse_interval(Real lower,Real upper = 2) const;
  FractionMeasure measure() const { return measure_; }
  Real cell_lower(std::size_t index) const;
  Real cell_upper(std::size_t index) const;
  const std::vector<Real>& fractions() const { return z_; }
  const std::vector<Real>& weights() const { return weights_; }
  Real minimum_z() const { return minimum_z_; }
private:
  friend class FractionPairWorkspace;
  Real minimum_z_,cell_fraction_,spacing_;
  FractionMeasure measure_;
  bool resolved_ = false;
  std::vector<Real> z_,weights_,prefix_,inverse_prefix_,upper_,lower_;
  Real tail(Real lower,bool inverse) const;
  std::array<Real,2> tails(Real lower) const;
  Real integrated_tail(Real product,Real first_lower,Real first_upper,Real second_lower) const;
};

// Per-evaluation scratch storage. Reuses exact cut/search stencils across
// signed table channels on compatible grids. No physics tables or global
// mutable cache; use a separate workspace on each concurrent worker.
class FractionPairWorkspace {
public:
  Real apply(const SignedFractionTable&,const SignedFractionTable&,
             Real lower1,Real lower2,Real product_min,Real product_max);
  void clear() { plans_.clear(); }
private:
  struct Query {std::size_t index=0;Real partial=0,inverse=0;bool zero=true;};
  struct Interval {Query lower,upper;std::size_t begin=0,end=0;bool atom=false,exclude_endpoint=false;};
  struct Integral {Query lower,upper;Real a=0,b=0,product=0;bool zero=true;};
  struct Row {std::size_t index;Interval interval;Integral lower,upper;Real density=0;};
  struct Plan {
    Real h1,h2,min1,min2,l1,l2,lo,hi;
    std::size_t n1,n2;
    FractionMeasure m1,m2;
    std::vector<Row> rows;
  };
  std::vector<Plan> plans_;
  static Query query(const SignedFractionTable&,Real);
  static std::array<Real,2> value(const SignedFractionTable&,const Query&);
  static Interval interval(const SignedFractionTable&,Real,Real);
  static Real value(const SignedFractionTable&,const Interval&);
  static Integral integral(const SignedFractionTable&,Real,Real,Real,Real);
  static Real value(const SignedFractionTable&,const Integral&);
};

// Born-parent FF convolution only: preserves host spin/recoil through its input
// momenta. Restores the continuum FF O(alpha) coefficient, NOT the missing
// process-specific finite hard coefficient. The two-leg Z action counts
// inclusive identified opposite-charge pairs from distinct hard directions.
class LabFragmentation {
public:
  explicit LabFragmentation(const QEDXSpaceEvolution&,Real minimum_z = 1e-4L,
                            FractionMeasure = FractionMeasure::atoms);
  // Fixed-alpha photonic U_LL only, with the same continuum D1 restoration.
  // Pairs, running and NLO kernels are absent: additive NLO+LL, not NLL.
  LabFragmentation(const LLNonSinglet&,const QEDParameters&,Real minimum_z = 1e-4L);
  // U + a U*d1, at most one boundary insertion globally in a two-leg action.
  explicit LabFragmentation(const LLBoundaryNonSinglet&,Real minimum_z = 1e-4L);
  // Fixed-alpha unitary photonic NLL, one boundary OR P1 insertion globally.
  explicit LabFragmentation(const PhotonicNLL&,Real minimum_z = 1e-4L);
  LabFFResult w(const FourVector& charged,const FourVector& neutrino,
                const WLabCuts& cuts = {}) const;
  LabFFResult z(const FourVector& first,const FourVector& second,
                const FiducialCuts& cuts = {}) const;
  const SignedFractionTable& same() const { return same_; }
  const SignedFractionTable& opposite() const { return opposite_; }
  Real a_initial() const { return a_; }
private:
  Real a_,mass_log_;
  SignedFractionTable same_,opposite_,c1_,c2_;
  std::optional<SignedFractionTable> boundary_,initial_boundary_;
  Real continuum_interval(Real lower,Real upper = 2) const;
  LabFFResult finish(LabFFResult) const;
};
} // namespace leptonff
