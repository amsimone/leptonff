#include "leptonff/xspace_evolution.hpp"
#include "fft.hpp"
#include "matrix_evolution.hpp"
#include "flavour_evolution.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <utility>

namespace leptonff {
namespace {
using Complex = detail::Complex;
using Spectrum = std::array<std::array<std::vector<Complex>,3>,3>;
struct Cells {
  std::vector<Real> weights;
  Real absolute_tail = 0;
};
Cells cells(const Distribution& d, const QEDXSpaceGrid& grid,
            const IntegrationOptions& options) {
  const auto n = grid.bins;
  const Real h = grid.spacing(), z1 = std::exp(-h), cap = std::exp(-h*(n-1));
  Cells result{std::vector<Real>(n,0),0};
  // First-cell test vanishes at the endpoint, so the full plus distribution
  // acts without a singular virtual term. Deposit its finite neighbour part.
  result.weights[1] = d.apply([=](Real z) { return z*(1-z)/(1-z1); },z1,options);
  for (std::size_t j = 1; j+1 < n; ++j) {
    const Real za = std::exp(-h*j), zb = std::exp(-h*(j+1));
    const Real width = za-zb;
    const auto density = [&](Real u) {
      const Real z = std::exp(-u);
      return z*z*d.resolved_density(z);
    };
    const Real total = integrate(density,j*h,(j+1)*h,options).value;
    const Real right = integrate([&](Real u) {
      const Real z = std::exp(-u);
      return density(u)*(za-z)/width;
    },j*h,(j+1)*h,options).value;
    result.weights[j] += total-right;
    result.weights[j+1] += right;
  }
  auto tail_options = options;
  tail_options.absolute_tolerance *= cap;
  const auto tail = integrate([&](Real t) {
    const Real z = cap*t;
    return cap*z*d.resolved_density(z);
  },0,1,tail_options);
  const auto absolute_tail = integrate([&](Real t) {
    const Real z = cap*t;
    return cap*z*std::abs(d.resolved_density(z));
  },0,1,tail_options);
  result.absolute_tail = absolute_tail.value+absolute_tail.estimated_error;
  Real resolved = 0;
  for (std::size_t j = 1; j < n; ++j) resolved += result.weights[j];
  // This is an individual distribution identity, not a row normalisation.
  // Subtract the omitted real tail explicitly instead of moving it to z=1.
  result.weights[0] = d.moment(2,options)-resolved-tail.value;
  return result;
}

Real row_norm(const MomentMatrix& a) {
  Real result = 0;
  for (const auto& row : a) result = std::max(result,row[0]+row[1]+row[2]);
  return result;
}
Real logarithmic_norm(const LatticeMatrix& kernel, Real h, Real exponent) {
  Real result = -1e300L;
  for (unsigned i = 0; i < 3; ++i) {
    Real row = kernel[i][i][0]; // retain signed virtual diagonal
    for (unsigned j = 0; j < 3; ++j) {
      if (i != j) row += std::abs(kernel[i][j][0]);
      for (std::size_t k = 1; k < kernel[i][j].size(); ++k)
        row += std::abs(kernel[i][j][k])*std::exp(exponent*h*k);
    }
    result = std::max(result,row);
  }
  return result;
}
Real boundary_norm(const LatticeMatrix& boundary, Real a, Real h, Real exponent) {
  MomentMatrix variations{};
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      variations[i][j] = std::abs((i==j ? 1 : 0)+a*boundary[i][j][0]);
      for (std::size_t k = 1; k < boundary[i][j].size(); ++k)
        variations[i][j] += a*std::abs(boundary[i][j][k])*std::exp(exponent*h*k);
    }
  return row_norm(variations);
}
Spectrum transform(const LatticeMatrix& a) {
  Spectrum spectrum;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      spectrum[i][j].assign(a[i][j].begin(),a[i][j].end());
      detail::fft(spectrum[i][j],false);
    }
  return spectrum;
}
detail::Matrix3<Complex> at(const Spectrum& a, std::size_t k) {
  detail::Matrix3<Complex> result;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) result[i][j] = a[i][j][k];
  return result;
}
void put(Spectrum& a, std::size_t k, const detail::Matrix3<Complex>& matrix,
         std::size_t n) {
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      a[i][j][k] = matrix[i][j];
      if (k && k != n/2) a[i][j][n-k] = std::conj(matrix[i][j]);
    }
}
Real inverse(Spectrum& a, LatticeMatrix& result) {
  Real imaginary = 0;
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) {
      detail::fft(a[i][j],true);
      result[i][j].resize(a[i][j].size());
      for (std::size_t k = 0; k < a[i][j].size(); ++k) {
        result[i][j][k] = a[i][j][k].real();
        imaginary += std::abs(a[i][j][k].imag());
      }
    }
  return imaginary;
}
} // namespace

void QEDXSpaceGrid::validate() const {
  if (bins < 16 || bins > (std::size_t{1}<<18) || (bins & (bins-1)))
    throw std::invalid_argument("Coupled grid bins must be a power of two in [16,2^18]");
  if (!std::isfinite(log_span) || log_span <= 0 || log_span > 80 ||
      !std::isfinite(wrap_exponent) || wrap_exponent <= 0 || wrap_exponent >= 1 ||
      !std::isfinite(max_wrap_bound) || max_wrap_bound <= 0 || max_wrap_bound >= 1 ||
      !std::isfinite(max_tail_bound) || max_tail_bound <= 0 || max_tail_bound >= 1)
    throw std::invalid_argument("Invalid coupled grid span/bound controls");
}

QEDXSpaceEvolution::QEDXSpaceEvolution(const MomentEvolutionParameters& parameters,
    const QEDXSpaceGrid& grid, const MomentEvolutionOptions& options)
    : parameters_(parameters),grid_(grid) {
  parameters.validate(); grid.validate(); options.validate();
  const Real L = 2*(std::log(parameters.mu_frag)-std::log(parameters.mu_0));
  if (L < 0) throw std::invalid_argument("Coupled x-space reference requires forward evolution");
  if(parameters.flavour_scheme) {
    auto r=detail::evolve_flavour_xspace(parameters,grid,options);
    weights_=std::move(r.weights);propagator_=std::move(r.propagator);kernel0_=std::move(r.kernel0);kernel1_=std::move(r.kernel1);
    boundary1_=std::move(r.boundary1);coefficient1_=std::move(r.c1);coefficient2_=std::move(r.c2);
    diagnostics_=r.diagnostics;flavour_weights_=std::move(r.full_weights);return;
  }
  const auto n = grid.bins;
  const Real h = grid.spacing(), a = parameters.alpha_msbar_at_mu0/(2*std::acos(-1.0L));
  const QEDParameters q{parameters.alpha_msbar_at_mu0,parameters.muon_mass,
                       parameters.mu_0,parameters.mu_frag};
  MomentMatrix tail0{},tail1{},tailb{};
  for (const auto& pair : {std::pair<unsigned,unsigned>{0,0},{0,1},{0,2},{2,0},{2,2}}) {
    const auto i = pair.first, j = pair.second;
    const auto p = one_family_species[i], o = one_family_species[j];
    auto make = [&](const Distribution& d,const char* label) {
      try { return cells(d,grid,options.kernel_integration); }
      catch(const std::exception& error) {
        throw std::runtime_error(std::string(label)+" "+species_name(p)+"->"+species_name(o)+": "+error.what());
      }
    };
    auto c0 = make(leading_kernel(p,o),"P0");
    auto c1 = make(next_to_leading_kernel(p,o),"P1");
    auto cb = make(initial_condition_coefficient(p,o,q),"d1");
    kernel0_[i][j]=std::move(c0.weights); kernel1_[i][j]=std::move(c1.weights);
    boundary1_[i][j]=std::move(cb.weights);
    tail0[i][j]=c0.absolute_tail; tail1[i][j]=c1.absolute_tail; tailb[i][j]=cb.absolute_tail;
  }
  for (auto* lattice : {&kernel0_,&kernel1_,&boundary1_}) {
    (*lattice)[1][1]=(*lattice)[0][0]; (*lattice)[1][0]=(*lattice)[0][1];
    (*lattice)[1][2]=(*lattice)[0][2]; (*lattice)[2][1]=(*lattice)[2][0];
  }
  for (auto* tail : {&tail0,&tail1,&tailb}) {
    (*tail)[1][1]=(*tail)[0][0]; (*tail)[1][0]=(*tail)[0][1];
    (*tail)[1][2]=(*tail)[0][2]; (*tail)[2][1]=(*tail)[2][0];
  }
  // The scalar running integrals are independent of the generator.
  std::size_t steps; Real shift;
  const auto running = detail::evolve_generators(MomentMatrix{},MomentMatrix{},a,L,
      parameters.coupling,options,steps,shift);
  const Real sigma = grid.wrap_exponent;
  diagnostics_.weighted_wrap_bound = std::exp(-sigma*grid.log_span
      +running.integral_a*logarithmic_norm(kernel0_,h,sigma)
      +running.integral_a2*logarithmic_norm(kernel1_,h,sigma))
      *boundary_norm(boundary1_,a,h,sigma);
  const Real growth = std::exp(running.integral_a*logarithmic_norm(kernel0_,h,0)
      +running.integral_a2*logarithmic_norm(kernel1_,h,0));
  const Real missing = running.integral_a*row_norm(tail0)+running.integral_a2*row_norm(tail1);
  diagnostics_.weighted_tail_bound = growth*(boundary_norm(boundary1_,a,h,0)*std::expm1(missing)
      +a*row_norm(tailb)*std::exp(missing));
  if (!std::isfinite(diagnostics_.weighted_wrap_bound) ||
      diagnostics_.weighted_wrap_bound > grid.max_wrap_bound ||
      !std::isfinite(diagnostics_.weighted_tail_bound) ||
      diagnostics_.weighted_tail_bound > grid.max_tail_bound)
    throw std::invalid_argument("Coupled log span fails weighted wrap/tail bounds");

  const auto f0 = transform(kernel0_), f1 = transform(kernel1_), fb = transform(boundary1_);
  Spectrum fe=f0,fd=f0,fc1=f0,fc2=f0;
  const Real beta0 = parameters.coupling == QEDCouplingMode::fixed ? 0 : 2.0L/3;
  for (std::size_t k = 0; k <= n/2; ++k) {
    const auto p0=at(f0,k),p1=at(f1,k),b1=at(fb,k);
    const auto state = detail::evolve_generators(p0,p1,a,L,parameters.coupling,
        options,steps,shift);
    diagnostics_.maximum_steps=std::max(diagnostics_.maximum_steps,steps);
    diagnostics_.maximum_frequency_step_shift=std::max(diagnostics_.maximum_frequency_step_shift,shift);
    auto boundary=detail::identity_matrix<Complex>();
    auto c1=b1,c2=p1;
    const auto square=detail::matrix_product(p0,p0),pb=detail::matrix_product(p0,b1);
    for (unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) {
      boundary[i][j]+=a*b1[i][j];
      c1[i][j]+=L*p0[i][j];
      c2[i][j]=L*p1[i][j]+L*pb[i][j]+L*L/2*(square[i][j]+beta0*p0[i][j]);
    }
    put(fe,k,state.e,n); put(fd,k,detail::matrix_product(state.e,boundary),n);
    put(fc1,k,c1,n); put(fc2,k,c2,n);
  }
  diagnostics_.weighted_step_shift_bound=3*std::sqrt(Real(n))
      *diagnostics_.maximum_frequency_step_shift*boundary_norm(boundary1_,a,h,0);
  diagnostics_.imaginary_mass=inverse(fd,weights_);
  inverse(fe,propagator_); inverse(fc1,coefficient1_); inverse(fc2,coefficient2_);
  for (const auto& row : weights_) for(const auto& channel:row) for(Real w:channel) {
    if (!std::isfinite(w)) throw std::runtime_error("Non-finite coupled lattice weight");
    if (w<0) diagnostics_.negative_weight_mass-=w;
  }
  if (!std::isfinite(diagnostics_.imaginary_mass) || diagnostics_.imaginary_mass>1e-8L)
    throw std::runtime_error("Coupled inverse FFT imaginary residual exceeds tolerance");
}

Real QEDXSpaceEvolution::z(std::size_t index) const {
  if (index>=grid_.bins) throw std::out_of_range("Coupled lattice index");
  return std::exp(-grid_.spacing()*index);
}
MomentMatrix QEDXSpaceEvolution::moment(Real N) const {
  if (!std::isfinite(N) || N<2)
    throw std::invalid_argument("Energy-weighted lattice moments require N>=2");
  MomentMatrix result{};
  for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j)
    for(std::size_t k=0;k<grid_.bins;++k)
      result[i][j]+=weights_[i][j][k]*std::exp(-(N-2)*grid_.spacing()*k);
  return result;
}
std::vector<std::array<Real,3>> QEDXSpaceEvolution::flavour_moment(Real N) const {
  if(!parameters_.flavour_scheme) {
    const auto m=moment(N);return {m[0],m[1],m[2]};
  }
  if(!std::isfinite(N) || N<2)throw std::invalid_argument("Energy-weighted flavour moments require N>=2");
  std::vector<std::array<Real,3>> result(flavour_weights_.size());
  for(std::size_t i=0;i<result.size();++i)for(unsigned j=0;j<3;++j)
    for(std::size_t k=0;k<grid_.bins;++k)result[i][j]+=flavour_weights_[i][j][k]*std::exp(-(N-2)*grid_.spacing()*k);
  return result;
}
Real QEDXSpaceEvolution::apply(Species parent,Species observed,const Function& measurement,
                              Real minimum_z) const {
  species_name(parent); species_name(observed);
  if (!std::isfinite(minimum_z) || minimum_z<=0 || minimum_z>1 ||
      minimum_z<z(grid_.bins-1))
    throw std::invalid_argument("Truncated action needs a positive lower bound inside the grid");
  const auto& weights=weights_[static_cast<unsigned>(parent)][static_cast<unsigned>(observed)];
  Real result=0;
  for(std::size_t k=0;k<grid_.bins;++k) {
    const Real fraction=z(k);
    if(fraction<minimum_z) break;
    const Real value=measurement(fraction);
    if(!std::isfinite(value)) throw std::invalid_argument("Non-finite coupled measurement");
    result+=weights[k]*value/fraction;
  }
  if(!std::isfinite(result)) throw std::runtime_error("Non-finite coupled action");
  return result;
}

} // namespace leptonff
