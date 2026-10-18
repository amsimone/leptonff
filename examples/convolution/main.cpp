#include <leptonff/core.hpp>
#include <iomanip>
#include <iostream>

int main() {
  using namespace leptonff;
  const auto ff=fragmentation_order1(Species::muon_minus,Species::muon_minus,{});
  const Function coefficient=[](Real x) { return x*(1-x); };
  std::cout << std::setprecision(12);
  for (Real x : {.1L,.5L,.9L})
    std::cout << "x=" << x << " convolution " << ff.convolve(coefficient,x) << '\n';
}
