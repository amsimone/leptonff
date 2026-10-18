#include "fft.hpp"
#include <algorithm>
#include <cmath>

namespace leptonff { namespace detail {
void fft(std::vector<Complex>& values, bool inverse) {
  const std::size_t n = values.size();
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(values[i], values[j]);
  }
  for (std::size_t length = 2; length <= n; length <<= 1) {
    const Real angle = (inverse ? 2 : -2)*std::acos(-1.0L)/length;
    const Complex step = std::polar(Real{1}, angle);
    for (std::size_t start = 0; start < n; start += length) {
      Complex phase = 1;
      for (std::size_t j = 0; j < length/2; ++j) {
        // Bound accumulated twiddle error on 53-bit long-double platforms.
        if (j % 256 == 0) phase = std::polar(Real{1}, angle*j);
        const Complex even = values[start+j];
        const Complex odd = phase*values[start+j+length/2];
        values[start+j] = even+odd;
        values[start+j+length/2] = even-odd;
        phase *= step;
      }
    }
  }
  if (inverse) for (auto& value : values) value /= n;
}
}} // namespace leptonff::detail
