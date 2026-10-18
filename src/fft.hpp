#pragma once
#include "leptonff/quadrature.hpp"
#include <complex>

namespace leptonff { namespace detail {
using Complex = std::complex<Real>;
void fft(std::vector<Complex>& values, bool inverse);
}} // namespace leptonff::detail
