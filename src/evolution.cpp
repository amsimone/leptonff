#include "leptonff/evolution.hpp"
#include "fft.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

namespace leptonff {
namespace {
using Complex = std::complex<Real>;
constexpr Real euler_gamma = .577215664901532860606512090082402431L;

struct Sum {
  Real value = 0, correction = 0;
  void add(Real x) {
    const Real corrected = x - correction;
    const Real next = value + corrected;
    correction = (next - value) - corrected;
    value = next;
  }
};

Real digamma(Real x) {
  Real shift = 0;
  while (x < 16) { shift -= 1/x; x += 1; }
  const Real r = 1/(x*x);
  return shift + std::log(x) - 1/(2*x)
      - r*(1.0L/12 - r*(1.0L/120 - r*(1.0L/252
      - r*(1.0L/240 - r*(1.0L/132)))));
}

// Integral of nu(u)=exp(-u)*(1+exp(-2u))/(1-exp(-u)). Evaluating
// expm1/log1p combinations avoids cancellation for both tiny and large u.
Real jump_integral(Real lower, Real upper) {
  const Real e = std::exp(-lower);
  const Real change = std::expm1(-(upper-lower));
  return 2*std::log1p(-e*change/(-std::expm1(-lower)))
      + e*change + e*e*std::expm1(-2*(upper-lower))/2;
}

Real jump_tail(Real lower) {
  const Real e = std::exp(-lower);
  return -2*std::log1p(-e) - e - e*e/2;
}

// Integral of the boundary measure nu(u)*(L0-1-2*log(1-exp(-u))).
// t=1-exp(-u); integrate (2/t-2+t)*log(t) with stable differences.
Real boundary_integral(Real lower, Real upper, Real L0) {
  const Real e = std::exp(-lower);
  const Real t0 = -std::expm1(-lower);
  const Real delta = -e*std::expm1(-(upper-lower));
  const Real t1 = t0+delta;
  const Real log0 = lower < .5L ? std::log(t0) : std::log1p(-e);
  const Real dlog = std::log1p(delta/t0);
  const Real log_integral = dlog*(2*log0+dlog)
      -2*(delta*(log0-1)+t1*dlog)
      +delta*(t1+t0)*(log0-.5L)/2+t1*t1*dlog/2;
  return (L0-1)*jump_integral(lower,upper)-2*log_integral;
}

Real lattice_action(const std::vector<Real>& weights, Real spacing,
                    const Function& measurement) {
  Sum result;
  for (std::size_t j = 0; j < weights.size(); ++j) {
    if (weights[j] == 0) continue;
    const Real value = measurement(std::exp(-spacing*j));
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite lattice measurement");
    result.add(weights[j]*value);
  }
  if (!std::isfinite(result.value)) throw std::runtime_error("Non-finite lattice action");
  return result.value;
}

using detail::fft;

} // namespace

LatticeConvolution circular_lattice_convolution(const std::vector<Real>& first,
                                               const std::vector<Real>& second) {
  const auto n = first.size();
  if (n < 16 || n > (std::size_t{1} << 22) || (n & (n-1)) != 0 || second.size() != n)
    throw std::invalid_argument("Incompatible periodic convolution arrays");
  std::vector<Complex> f(n), g(n);
  Sum sum_f, sum_g, variation_f, variation_g;
  for (std::size_t j = 0; j < n; ++j) {
    if (!std::isfinite(first[j]) || !std::isfinite(second[j]))
      throw std::invalid_argument("Non-finite periodic convolution weights");
    f[j] = first[j];
    g[j] = second[j];
    sum_f.add(first[j]); sum_g.add(second[j]);
    variation_f.add(std::abs(first[j])); variation_g.add(std::abs(second[j]));
  }
  fft(f,false); fft(g,false);
  for (std::size_t j = 0; j < n; ++j) f[j] *= g[j];
  fft(f,true);
  LatticeConvolution result;
  result.weights.resize(n);
  Sum number, imaginary;
  for (std::size_t j = 0; j < n; ++j) {
    result.weights[j] = f[j].real();
    number.add(result.weights[j]);
    imaginary.add(std::abs(f[j].imag()));
  }
  result.number_residual = number.value-sum_f.value*sum_g.value;
  result.imaginary_weight_mass = imaginary.value;
  const Real tolerance = 1e-9L*std::max(Real{1},variation_f.value*variation_g.value);
  if (!std::isfinite(result.number_residual) ||
      std::abs(result.number_residual) > tolerance || imaginary.value > tolerance)
    throw std::runtime_error("Periodic convolution FFT residual exceeds tolerance");
  return result;
}

void LLGrid::validate() const {
  if (bins < 16 || bins > (std::size_t{1} << 22) || (bins & (bins-1)) != 0)
    throw std::invalid_argument("LL bins must be a power of two in [16,2^22]");
  if (!std::isfinite(log_span) || log_span <= 0 || log_span > 600)
    throw std::invalid_argument("LL log_span must be finite and in (0,600]");
  if (!std::isfinite(soft_log_cutoff) || soft_log_cutoff < spacing()/2 ||
      soft_log_cutoff >= log_span-spacing()/2)
    throw std::invalid_argument("LL soft cutoff must be >= half a grid cell and inside the grid");
  if (!std::isfinite(max_wrap_probability) || max_wrap_probability <= 0 ||
      max_wrap_probability >= 1)
    throw std::invalid_argument("Invalid LL wrap-probability tolerance");
}

Real LLGrid::spacing() const { return log_span/bins; }

Real non_singlet_kernel_moment(Real N) {
  if (!std::isfinite(N) || N <= 0)
    throw std::invalid_argument("Non-singlet kernel moment requires finite N>0");
  if (N == 1) return 0;
  return 1.5L - 2*(digamma(N+1)+euler_gamma) + 1/(N*(N+1));
}

LLNonSinglet::LLNonSinglet(Real time, const LLGrid& grid)
    : time_(time), grid_(grid) {
  grid_.validate();
  if (!std::isfinite(time_) || time_ < 0)
    throw std::invalid_argument("LL evolution time must be finite and nonnegative");
  weights_.resize(grid_.bins, 0);

  const Real h = grid_.spacing();
  std::vector<Complex> frequencies(grid_.bins, 0);
  Sum rate, tilted_rate;
  constexpr Real tilt = .9L;
  for (std::size_t j = 1; j < grid_.bins; ++j) {
    const Real lower = std::max(grid_.soft_log_cutoff, (j-.5L)*h);
    const Real upper = (j+.5L)*h;
    if (lower >= upper) continue;
    const Real weight = jump_integral(lower, upper);
    if (!std::isfinite(weight) || weight < 0)
      throw std::runtime_error("Invalid integrated LL jump rate");
    frequencies[j] = weight;
    rate.add(weight);
    tilted_rate.add(weight*std::expm1(tilt*j*h));
  }
  diagnostics_.resolved_jump_rate = rate.value;
  diagnostics_.no_resolved_emission_weight = std::exp(-time_*rate.value);
  diagnostics_.wrap_probability_bound = std::exp(time_*tilted_rate.value-tilt*grid_.log_span);
  diagnostics_.omitted_large_jump_probability_bound = time_*jump_tail((grid_.bins-.5L)*h);
  diagnostics_.mean_log_loss_discretisation_bound =
      2*time_*grid_.soft_log_cutoff + time_*rate.value*h/2;
  if (!std::isfinite(diagnostics_.wrap_probability_bound) ||
      diagnostics_.wrap_probability_bound > grid_.max_wrap_probability)
    throw std::invalid_argument("LL log_span too short: periodic-wrap probability exceeds tolerance");

  // Negative virtual rate at the endpoint implements the whole-plus kernel.
  frequencies[0] = -rate.value;
  generator_.resize(grid_.bins);
  for(std::size_t j=0;j<grid_.bins;++j)generator_[j]=frequencies[j].real();
  if (time_ == 0) { weights_[0] = 1; return; }
  fft(frequencies, false);
  for (auto& value : frequencies) value = std::exp(time_*value);
  frequencies[0] = 1; // exact number-conservation eigenvalue
  fft(frequencies, true);
  Sum normalisation, negatives, imaginary;
  for (std::size_t j = 0; j < grid_.bins; ++j) {
    weights_[j] = frequencies[j].real();
    normalisation.add(weights_[j]);
    if (weights_[j] < 0) negatives.add(-weights_[j]);
    imaginary.add(std::abs(frequencies[j].imag()));
  }
  diagnostics_.normalisation_residual = normalisation.value-1;
  diagnostics_.negative_weight_mass = negatives.value;
  diagnostics_.imaginary_weight_mass = imaginary.value;
  // Retain signed roundoff; do not hide it by clipping or renormalising.
  if (std::abs(diagnostics_.normalisation_residual) > 1e-9L ||
      diagnostics_.negative_weight_mass > 1e-9L ||
      diagnostics_.imaginary_weight_mass > 1e-9L)
    throw std::runtime_error("LL FFT numerical residuals exceed tolerance");
}

Real LLNonSinglet::z(std::size_t index) const {
  if (index >= weights_.size()) throw std::out_of_range("LL grid index");
  return std::exp(-grid_.spacing()*index);
}

Real LLNonSinglet::apply(const Function& measurement) const {
  Sum result;
  for (std::size_t j = 0; j < weights_.size(); ++j) {
    if (weights_[j] == 0) continue;
    const Real value = measurement(z(j));
    if (!std::isfinite(value)) throw std::invalid_argument("Non-finite LL measurement");
    result.add(weights_[j]*value);
  }
  if (!std::isfinite(result.value)) throw std::runtime_error("Non-finite LL action");
  return result.value;
}

Real LLNonSinglet::moment(Real N) const {
  if (!std::isfinite(N) || N < 1)
    throw std::invalid_argument("Numerical LL moments require finite N>=1");
  return apply([N](Real fraction) { return std::pow(fraction, N-1); });
}

Real LLNonSinglet::exact_moment(Real N) const {
  const Real value = std::exp(time_*non_singlet_kernel_moment(N));
  if (!std::isfinite(value)) throw std::overflow_error("Continuum LL moment overflow");
  return value;
}

LLBoundaryNonSinglet::LLBoundaryNonSinglet(const QEDParameters& parameters,
                                          const LLGrid& grid)
    : parameters_(parameters),
      evolution_(parameters.expansion_parameter()
          *2*(std::log(parameters.mu_frag)-std::log(parameters.mu_0)),grid) {
  const Real L0 = 2*(std::log(parameters_.mu_0)-std::log(parameters_.muon_mass));
  const Real L = 2*(std::log(parameters_.mu_frag)-std::log(parameters_.mu_0));
  const Real h = grid.spacing();
  initial_.resize(grid.bins,0);
  order1_.resize(grid.bins,0);
  Sum initial_rate, order1_rate, tilted_variation;
  for (std::size_t j = 1; j < grid.bins; ++j) {
    const Real lower = std::max(grid.soft_log_cutoff,(j-.5L)*h);
    const Real upper = (j+.5L)*h;
    if (lower >= upper) continue;
    initial_[j] = boundary_integral(lower,upper,L0);
    order1_[j] = initial_[j]+L*jump_integral(lower,upper);
    initial_rate.add(initial_[j]);
    order1_rate.add(order1_[j]);
    tilted_variation.add(std::abs(initial_[j])*std::exp(.9L*j*h));
  }
  initial_[0] = -initial_rate.value;
  order1_[0] = -order1_rate.value;
  tilted_variation.add(std::abs(initial_[0]));
  diagnostics_.wrap_coefficient_bound = evolution_.diagnostics().wrap_probability_bound
      *tilted_variation.value;
  const Real end = (grid.bins-.5L)*h;
  diagnostics_.omitted_tail_coefficient_bound = jump_tail(end)
      *(std::abs(L0-1)-2*std::log1p(-std::exp(-end)));
  if (parameters_.expansion_parameter()*diagnostics_.wrap_coefficient_bound >
      grid.max_wrap_probability)
    throw std::invalid_argument("LL boundary wrap bound exceeds tolerance: increase log_span");

  std::vector<Complex> boundary(grid.bins), evolution(grid.bins);
  for (std::size_t j = 0; j < grid.bins; ++j) {
    boundary[j] = initial_[j];
    evolution[j] = evolution_.weights()[j];
  }
  fft(boundary,false);
  fft(evolution,false);
  for (std::size_t j = 0; j < grid.bins; ++j) boundary[j] *= evolution[j];
  boundary[0] = 0; // exact boundary number-conservation eigenvalue
  fft(boundary,true);
  evolved_.resize(grid.bins);
  Sum initial_sum, evolved_sum, imaginary;
  for (std::size_t j = 0; j < grid.bins; ++j) {
    evolved_[j] = boundary[j].real();
    initial_sum.add(initial_[j]);
    evolved_sum.add(evolved_[j]);
    imaginary.add(std::abs(boundary[j].imag()));
  }
  diagnostics_.initial_number_residual = initial_sum.value;
  diagnostics_.evolved_number_residual = evolved_sum.value;
  diagnostics_.imaginary_weight_mass = imaginary.value;
  if (std::abs(initial_sum.value) > 1e-8L || std::abs(evolved_sum.value) > 1e-8L ||
      imaginary.value > 1e-8L)
    throw std::runtime_error("Boundary coefficient FFT residuals exceed tolerance");
}

Real LLBoundaryNonSinglet::apply(const Function& measurement) const {
  return evolution_.apply(measurement)+parameters_.expansion_parameter()
      *lattice_action(evolved_,evolution_.grid().spacing(),measurement);
}

Real LLBoundaryNonSinglet::boundary_moment(Real N) const {
  if (!std::isfinite(N) || N < 1)
    throw std::invalid_argument("Numerical LL boundary moments require finite N>=1");
  return lattice_action(evolved_,evolution_.grid().spacing(),
      [N](Real z) { return std::pow(z,N-1); });
}

Real LLBoundaryNonSinglet::exact_boundary_moment(Real N,
                                               const IntegrationOptions& options) const {
  return evolution_.exact_moment(N)*initial_condition_coefficient(
      Species::muon_minus,Species::muon_minus,parameters_).moment(N,options);
}

Real LLBoundaryNonSinglet::first_order_apply(const Function& measurement) const {
  return lattice_action(order1_,evolution_.grid().spacing(),measurement);
}
} // namespace leptonff
