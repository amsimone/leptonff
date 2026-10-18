#include <leptonff/core.hpp>
#include <iomanip>
#include <iostream>

int main() {
  using namespace leptonff;
  const QEDParameters input; // one-family MSbar, alpha/(2*pi) expansion
  const auto ff=fragmentation_order1(Species::muon_minus,Species::muon_minus,input);
  std::cout << std::setprecision(12)
            << "N=2 fragmentation moment " << ff.moment(2) << '\n'
            << "Endpoint bin [0.9,1] " << ff.integral(.9L,1) << '\n';
  EvolutionTheory theory;
  theory.alpha=AlphaInput::on_shell(input.alpha);
  theory.flavours=QEDFlavourScheme::leptons();
  const auto result=evolve_qed_moment(2,theory.parameters());
  std::cout << "Coupled NLL N=2 moment " << result.fragmentation[0][0] << '\n';
}
