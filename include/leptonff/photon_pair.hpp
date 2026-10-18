#pragma once

#include "leptonff/fiducial_vector.hpp"

namespace leptonff {
// One identified photon and ONE specified primary charge in V -> mu- mu+ gamma.
// y = 2*k_gamma.k_mu/M^2 (massless pair invariant / M^2). Coefficients of
// a=alpha/(2*pi), normalised to massless vector Born; no factor of two.
// MSbar: subtract Born x D1[mu -> gamma]. No virtual or endpoint delta.
// The unweighted small-y integral diverges; Mellin moments require N>1.
Distribution vector_photon_pair_mass_hard_coefficient(Real mass, Real mu_frag);

// Exact massive real spectrum in y=(s_gamma_mu-m_mu^2)/M^2, one primary charge.
// Returns zero beyond the physical endpoint 1-2*m_mu/M; requires 0<y<1.
Real massive_vector_photon_pair_spectrum(Real mass, Real muon_mass, Real y);

// y=s_gamma_mu/M^2, t=s_gamma_other_mu/M^2; y>0, t>=0, y+t<=1.
// Fractions rescale the photon and selected muon directions respectively.
// This measures a candidate pair AFTER photon conversion, not a photon veto.
Real vector_photon_pair_acceptance(Real mass, Real y, Real t,
    Real photon_fraction, Real muon_fraction, const FiducialCuts& cuts,
    const SingleLeptonSelection& selection = {});

struct PhotonPairHardCoefficient {
  Real projected = 0, geometry = 0, total = 0;
};

// H1[gamma,mu] for fixed fragmentation fractions, including local collinear
// subtraction. Positive pair-mass lower cut required: pairs in one unresolved
// direction are outside this two-distinct-parent measurement. The collinear
// control omits t/(1-y)<delta from an already finite geometry integral.
// soft_cutoff is not used: mass_min>0 excludes soft photons exactly.
PhotonPairHardCoefficient fiducial_vector_photon_pair_hard(Real mass,
    Real mu_frag, const FiducialCuts& cuts, Real photon_fraction = 1,
    Real muon_fraction = 1, const HardDecayOptions& options = {},
    const SingleLeptonSelection& selection = {});
} // namespace leptonff
