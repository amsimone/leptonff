#include <leptonff/core.hpp>
#include <iomanip>
#include <iostream>

int main() {
  using namespace leptonff;
  const auto ff=fragmentation_order1(Species::muon_minus,Species::muon_minus,{});
  const auto kernel=leading_kernel(Species::muon_minus,Species::muon_minus);
  std::cout << std::setprecision(18) << "{\n"
    << "  \"ff_N1\": " << ff.moment(1) << ",\n"
    << "  \"ff_N2\": " << ff.moment(2) << ",\n"
    << "  \"ff_N3\": " << ff.moment(3) << ",\n"
    << "  \"kernel_bin\": " << kernel.integral(.2L,.7L) << ",\n"
    << "  \"kernel_endpoint_bin\": " << kernel.integral(.7L,1) << ",\n"
    << "  \"kernel_convolution\": " << kernel.convolve([](Real) { return 1; },.3L) << "\n}\n";
}
