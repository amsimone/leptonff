#pragma once
#include "leptonff/production_decay.hpp"

namespace leptonff {
// Angular density with measure dOmega/(4*pi):
// 3/4 { U(1+c^2) + H0/2(1-3c^2) + H1 sin2theta cosphi
// + H2/2 sin^2theta cos2phi + H3 sintheta cosphi + H4 c
// + H5 sin^2theta sin2phi + H6 sin2theta sinphi + H7 sintheta sinphi }.
// H_i = U*A_i, kept UNNORMALISED so zero/negative matching terms are safe.
struct DrellYanAngularWeights {Real unpolarised=0;std::array<Real,8> angular{};};
// For neutral currents only the parity-even part is returned; this is valid
// for the charge-symmetric Z observables supported by DecayCalculator.
SpinTensor production_spin_tensor(SpinDecay,const DrellYanAngularWeights&,
                                 const SpatialRotation& angular_to_canonical);
// Explicit convention: z || unit(P1*)-unit(P2*), x || unit(P1*)+unit(P2*),
// y=z cross x. Columns map these axes into canonical LAB-boost rest axes.
// At qT=0, x is the projected LAB x direction. A backend must VERIFY this
// convention against its own recoil/azimuth convention before use.
SpatialRotation collins_soper_rotation(const FourVector& parent);
} // namespace leptonff
