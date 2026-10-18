#pragma once
#include "leptonff/lab_fragmentation.hpp"
#include <array>
#include <complex>

namespace leptonff {
enum class SpinDecay { wminus, wplus, neutral_vector };
using SpinTensor = std::array<std::array<std::complex<Real>,3>,3>;
using SpatialRotation = std::array<std::array<Real,3>,3>;

// Rest-frame Cartesian tensors. B is normalised to <B(n)>_sphere = I.
// Tr(rho B) is the angular density, with dOmega/(4*pi).
SpinTensor spin_born_tensor(SpinDecay,const std::array<Real,3>& direction);
// Signed dual angular projector: <Tr(rho B(n)) K(n)> = rho.
// Neutral-vector K reconstructs only the real symmetric part. Its use with
// a gamma*/Z host is restricted to charge-symmetric lepton measurements.
SpinTensor spin_angular_projector(SpinDecay,const std::array<Real,3>& direction);
Real spin_contract(const SpinTensor& density,const SpinTensor& decay);
Real spin_trace(const SpinTensor&);
SpinTensor rotate_spin_tensor(const SpinTensor&,const SpatialRotation&);
FourVector rotate_momentum(const FourVector&,const SpatialRotation&);
FourVector rest_frame_boost(const FourVector&,const FourVector& parent,int sign);
SpatialRotation decay_rotation(Real cosine,Real phi,Real twist);

struct SpinRadiativePoint {
  // W: charged lepton and neutrino. Neutral: negative and positive lepton.
  FourVector first,second,photon;
  SpinTensor real;
};
// Photon along +z; c is the first lepton cosine in the two-fermion rest
// frame. real is dC1/(dz dc dR), C1 multiplies alpha/(2*pi), dR is normalised
// Haar measure. z=(p_first+p_second)^2/M^2. Exact massive real current;
// W includes the charged-W emission graph (unitary-gauge propagator).
SpinRadiativePoint spin_radiative_point(SpinDecay,Real M,Real m,Real z,Real c);
Real spin_photon_ward_residual(SpinDecay,Real M,Real m,Real z,Real c);

// Unitarily normalised decay: inclusive C1=0. Virtual finite normalisation
// is transferred to the branching ratio. Kernels below are whole-plus,
// per charged leg; the W kernel includes its additional finite W term.
Real projected_decay_hard_kernel(SpinDecay,Real z,Real M,Real mu_frag);
Real w_projected_decay_hard(const FourVector&,const FourVector&,
                           Real mu_frag,const WLabCuts& = {});
Real z_projected_decay_hard(const FourVector&,const FourVector&,
                           Real mu_frag,const FiducialCuts& = {});

using DecayMeasurement = std::function<Real(const FourVector&,const FourVector&)>;
// Finite, leading-power angular remainder at a single integration point.
// Integrate dz dc dR. Only m=0 is used: neither a physical mass regulator nor
// a fitted generator correction enters H1. Remainder cutoffs must be refined.
SpinTensor decay_hard_geometry(SpinDecay,Real M,Real z,Real c,
    const SpatialRotation&,const FourVector& parent,const DecayMeasurement&,
    Real orientation_independent_reference = 0);
// Evaluate the bare operator and G[F_evolved-F_bare] with ONE shared current,
// rotation and projection. The difference vanishes pointwise at identity
// evolution. Evolved measurement must contain U only (no second finite term).
std::array<SpinTensor,2> decay_hard_geometry_pair(SpinDecay,Real M,Real z,Real c,
    const SpatialRotation&,const FourVector& parent,const DecayMeasurement& bare,
    const DecayMeasurement& evolved,Real orientation_independent_reference = 0);
} // namespace leptonff
