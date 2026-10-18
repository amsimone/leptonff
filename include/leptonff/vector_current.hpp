#pragma once

#include "leptonff/decay.hpp"

namespace leptonff {
struct VectorRadiativePoint {
  FourVector minus, plus, photon;
};

// Parent at rest, photon on +z; cos_theta is the mu- angle to +z in the
// dimuon rest frame. Global orientations are integrated in the density below.
// m=0 is supported away from cos_theta=+/-1; its angular integral diverges.
VectorRadiativePoint vector_radiative_point(Real parent_mass, Real muon_mass,
                                           Real z, Real cos_theta);

// Independent Dirac-trace real-emission matrix element, coefficient of
// alpha/(2*pi), relative to MASSLESS vector Born: dC_real/(dz dcos_theta).
Real vector_real_density(Real parent_mass, Real muon_mass, Real z, Real cos_theta);
Real vector_real_angular_integral(Real parent_mass, Real muon_mass, Real z,
                                 const IntegrationOptions& options = {});
// Max absolute matrix element after replacing the photon polarisation by k.
Real vector_current_ward_residual(Real parent_mass, Real muon_mass,
                                  Real z, Real cos_theta);
} // namespace leptonff
