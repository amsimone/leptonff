#pragma once
#include "leptonff/moment_evolution.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace leptonff { namespace detail {
template<class T> using Matrix3 = std::array<std::array<T,3>,3>;
template<class T> Matrix3<T> identity_matrix() {
  Matrix3<T> result{};
  for (unsigned i=0;i<3;++i) result[i][i]=1;
  return result;
}
template<class T> Matrix3<T> matrix_product(const Matrix3<T>& a,const Matrix3<T>& b) {
  Matrix3<T> result{};
  for (unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j)
    for(unsigned k=0;k<3;++k) result[i][j]+=a[i][k]*b[k][j];
  return result;
}
template<class Scalar> struct EvolutionState {
  Real a = 0, integral_a = 0, integral_a2 = 0;
  Matrix3<Scalar> e{};
};
template<class Scalar> EvolutionState<Scalar> add_scaled(EvolutionState<Scalar> y, const EvolutionState<Scalar>& k, Real scale) {
  y.a += scale*k.a;
  y.integral_a += scale*k.integral_a;
  y.integral_a2 += scale*k.integral_a2;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) y.e[i][j] += scale*k.e[i][j];
  return y;
}
template<class Scalar> EvolutionState<Scalar> derivative(const EvolutionState<Scalar>& y, const Matrix3<Scalar>& p0, const Matrix3<Scalar>& p1,
                 QEDCouplingMode mode) {
  if (!std::isfinite(y.a) || y.a < 0 || y.a >= 1/(2*std::acos(-1.0L)))
    throw std::runtime_error("Running alpha left the supported perturbative range [0,1)");
  EvolutionState<Scalar> d;
  if (mode != QEDCouplingMode::fixed) d.a = 2*y.a*y.a/3;
  if (mode == QEDCouplingMode::two_loop) d.a += y.a*y.a*y.a;
  d.integral_a = y.a;
  d.integral_a2 = y.a*y.a;
  Matrix3<Scalar> generator{};
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      generator[i][j] = y.a*p0[i][j]+y.a*y.a*p1[i][j];
  d.e = matrix_product(generator,y.e);
  return d;
}
template<class Scalar> EvolutionState<Scalar> integrate_steps(Real a, Real L, std::size_t steps, const Matrix3<Scalar>& p0,
                      const Matrix3<Scalar>& p1, QEDCouplingMode mode) {
  EvolutionState<Scalar> y{a,0,0,identity_matrix<Scalar>()};
  const Real h = L/steps;
  for (std::size_t step = 0; step < steps; ++step) {
    const EvolutionState<Scalar> k1 = derivative(y,p0,p1,mode);
    const EvolutionState<Scalar> k2 = derivative(add_scaled(y,k1,h/2),p0,p1,mode);
    const EvolutionState<Scalar> k3 = derivative(add_scaled(y,k2,h/2),p0,p1,mode);
    const EvolutionState<Scalar> k4 = derivative(add_scaled(y,k3,h),p0,p1,mode);
    y = add_scaled(add_scaled(add_scaled(add_scaled(y,k1,h/6),k2,h/3),k3,h/3),k4,h/6);
  }
  derivative(y,p0,p1,mode); // validate the final coupling as well
  return y;
}
template<class Scalar> Real state_distance(const EvolutionState<Scalar>& a, const EvolutionState<Scalar>& b) {
  Real difference = std::max({std::abs(a.a-b.a),std::abs(a.integral_a-b.integral_a),
                              std::abs(a.integral_a2-b.integral_a2)});
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      if (!std::isfinite(std::abs(a.e[i][j])) || !std::isfinite(std::abs(b.e[i][j])))
        throw std::runtime_error("Non-finite moment evolution matrix");
      difference = std::max(difference,std::abs(a.e[i][j]-b.e[i][j]));
    }
  return difference;
}

template<class T> EvolutionState<T> evolve_generators(const Matrix3<T>& p0,
    const Matrix3<T>& p1, Real a0, Real L, QEDCouplingMode mode,
    const MomentEvolutionOptions& options, std::size_t& steps, Real& shift) {
  EvolutionState<T> state{a0,0,0,identity_matrix<T>()};
  steps=0; shift=0;
  if (a0==0 || L==0) return state;
  steps=options.initial_steps;
  auto coarse=integrate_steps(a0,L,steps,p0,p1,mode);
  while (true) {
    if (steps>options.max_steps/2)
      throw std::runtime_error("Matrix evolution failed step refinement tolerance");
    steps*=2;
    state=integrate_steps(a0,L,steps,p0,p1,mode);
    shift=state_distance(state,coarse);
    if (shift<=options.step_tolerance) return state;
    coarse=state;
  }
}
}} // namespace leptonff::detail
