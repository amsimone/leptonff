#include "leptonff/vector_current.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace leptonff {
namespace {
using Complex = std::complex<Real>;
using Matrix = std::array<Complex,16>;

Matrix identity() {
  Matrix result{};
  for (int i = 0; i < 4; ++i) result[5*i] = 1;
  return result;
}
Matrix scale(Matrix m, Real factor) {
  for (auto& x : m) x *= factor;
  return m;
}
Matrix add(const Matrix& a, const Matrix& b) {
  Matrix result{};
  for (int i = 0; i < 16; ++i) result[i] = a[i]+b[i];
  return result;
}
Matrix product(const Matrix& a, const Matrix& b) {
  Matrix result{};
  for (int i = 0; i < 4; ++i)
    for (int k = 0; k < 4; ++k) {
      const auto x = a[4*i+k];
      if (x == Complex{}) continue;
      for (int j = 0; j < 4; ++j) result[4*i+j] += x*b[4*k+j];
    }
  return result;
}
Matrix adjoint(const Matrix& a) {
  Matrix result{};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) result[4*i+j] = std::conj(a[4*j+i]);
  return result;
}
const std::array<Matrix,4>& gamma() {
  static const auto matrices = [] {
    std::array<Matrix,4> g{};
    g[0][0] = g[0][5] = 1;
    g[0][10] = g[0][15] = -1;
    const std::array<std::array<Complex,4>,3> pauli{{
        {{0,1,1,0}},{{0,Complex{0,-1},Complex{0,1},0}},{{1,0,0,-1}}}};
    for (int a = 0; a < 3; ++a)
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) {
          g[a+1][4*i+j+2] = pauli[a][2*i+j];
          g[a+1][4*(i+2)+j] = -pauli[a][2*i+j];
        }
    return g;
  }();
  return matrices;
}
Matrix slash(const FourVector& p) {
  return add(add(scale(gamma()[0],p.energy),scale(gamma()[1],-p.px)),
             add(scale(gamma()[2],-p.py),scale(gamma()[3],-p.pz)));
}
Real dot(const FourVector& a, const FourVector& b) {
  return a.energy*b.energy-a.px*b.px-a.py*b.py-a.pz*b.pz;
}

Matrix current(const VectorRadiativePoint& p, const FourVector& epsilon,
               Real dminus, Real dplus, int mu) {
  // On-shell reduced FSR current, including emission from both charged legs.
  const auto first = add(scale(identity(),2*dot(p.minus,epsilon)),
                         product(slash(epsilon),slash(p.photon)));
  const auto second = add(scale(identity(),2*dot(p.plus,epsilon)),
                          product(slash(p.photon),slash(epsilon)));
  return add(scale(product(first,gamma()[mu]),1/dminus),
             scale(product(gamma()[mu],second),-1/dplus));
}

void validate(Real mass, Real muon, Real z, Real cosine) {
  if (!std::isfinite(mass) || !std::isfinite(muon) || mass <= 0 || muon < 0 ||
      muon >= mass/2 || !std::isfinite(z) || z <= 4*std::pow(muon/mass,2) ||
      z >= 1 || !std::isfinite(cosine) || std::abs(cosine) > 1 ||
      (muon == 0 && std::abs(cosine) == 1))
    throw std::invalid_argument("Invalid massive vector real-emission point");
}
} // namespace

VectorRadiativePoint vector_radiative_point(Real mass, Real muon, Real z, Real c) {
  validate(mass,muon,z,c);
  const Real root = std::sqrt(z);
  const Real energy = mass*root/2;
  const Real momentum = std::sqrt((energy-muon)*(energy+muon));
  const Real boost = (1+z)/(2*root), velocity = -(1-z)/(1+z);
  const Real px = momentum*std::sqrt((1-c)*(1+c));
  return {{boost*(energy+velocity*momentum*c),px,0,boost*(momentum*c+velocity*energy)},
          {boost*(energy-velocity*momentum*c),-px,0,boost*(-momentum*c+velocity*energy)},
          {mass*(1-z)/2,0,0,mass*(1-z)/2}};
}

Real vector_real_density(Real mass, Real muon, Real z, Real c) {
  const auto p = vector_radiative_point(mass,muon,z,c);
  const Real beta = std::sqrt(1-4*std::pow(muon/mass,2)/z);
  // Analytic invariants avoid subtracting nearly parallel lab-frame vectors.
  const Real dminus = mass*mass*(1-z)*(1-beta*c)/2;
  const Real dplus = mass*mass*(1-z)*(1+beta*c)/2;
  const auto minus = add(slash(p.minus),scale(identity(),muon));
  const auto plus = add(slash(p.plus),scale(identity(),-muon));
  Complex trace = 0;
  for (const FourVector epsilon : {FourVector{0,1,0,0},FourVector{0,0,1,0}})
    for (int mu = 1; mu <= 3; ++mu) {
      const auto amplitude = current(p,epsilon,dminus,dplus,mu);
      const auto conjugate = product(product(gamma()[0],adjoint(amplitude)),gamma()[0]);
      const auto square = product(product(product(minus,amplitude),plus),conjugate);
      for (int i = 0; i < 4; ++i) trace += square[5*i];
    }
  // Factorised three-body phase space / massless Born. The vector spin average
  // cancels between numerator and denominator; no photon averaging is made.
  const Real value = (1-z)*beta*trace.real()/16;
  if (!std::isfinite(value) || value < 0 ||
      std::abs(trace.imag()) > 1e-8L*std::max(Real{1},std::abs(trace.real())))
    throw std::runtime_error("Invalid vector-current trace");
  return value;
}

Real vector_real_angular_integral(Real mass, Real muon, Real z,
                                 const IntegrationOptions& options) {
  validate(mass,muon,z,0);
  if (muon == 0) throw std::invalid_argument("Massless real angular integral diverges");
  const Real ratio = 4*std::pow(muon/mass,2)/z;
  const Real beta = std::sqrt(1-ratio);
  const Real edge = std::log((1+beta)/std::sqrt(ratio)); // atanh(beta), stably
  // c=tanh(y)/beta resolves both mass-regulated collinear peaks.
  return integrate([&](Real y) {
    const Real cosine = std::clamp(std::tanh(y)/beta,Real{-1},Real{1});
    const Real jacobian = 1/(beta*std::pow(std::cosh(y),2));
    return vector_real_density(mass,muon,z,cosine)*jacobian;
  },-edge,edge,options,{0}).value;
}

Real vector_current_ward_residual(Real mass, Real muon, Real z, Real c) {
  const auto p = vector_radiative_point(mass,muon,z,c);
  // Use identical contractions in denominators for the algebraic Ward test.
  const Real dminus = 2*dot(p.minus,p.photon), dplus = 2*dot(p.plus,p.photon);
  Real norm = 0;
  for (int mu = 0; mu <= 3; ++mu)
    for (const auto element : current(p,p.photon,dminus,dplus,mu))
      norm = std::max(norm,std::abs(element));
  return norm;
}
} // namespace leptonff
