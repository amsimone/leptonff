#include "leptonff/core.hpp"
#include "check.hpp"

int main() {
  using namespace leptonff;
  EvolutionTheory theory;
  const auto legacy=MomentEvolutionParameters{};
  close(evolve_qed_moment(2,theory.parameters()).fragmentation[0][0],
        evolve_qed_moment(2,legacy).fragmentation[0][0],1e-14L,"explicit theory preserves backend");
  theory.flavours=QEDFlavourScheme::leptons();
  theory.alpha=AlphaInput::on_shell(1/137.035999084L);
  close(theory.parameters().alpha_msbar_at_mu0,
        leptonic_alpha_msbar(theory.alpha.value,theory.initial_scale,*theory.flavours),
        1e-16L,"explicit on-shell conversion");
  theory.alpha=AlphaInput::msbar(.0075L,2);
  throws<std::invalid_argument>([&]{ theory.parameters(); },"refuse silent reference-scale change");
  theory.initial_scale=2;
  close(theory.parameters().alpha_msbar_at_mu0,.0075L,1e-16L,"explicit MSbar reference");
  theory.alpha=AlphaInput::on_shell(.0073L);
  theory.flavours=QEDFlavourScheme::hadronic_pointlike();
  throws<std::invalid_argument>([&]{ theory.parameters(); },"no implicit hadronic OS conversion");
  AlphaInput bad{.007L,AlphaScheme::on_shell,1};
  throws<std::invalid_argument>([&]{ bad.validate(); },"on-shell scale ambiguity");

  const auto p=leading_kernel(Species::muon_minus,Species::muon_minus);
  const Real x=.3L;
  close(p.convolve([](Real){return 1;},x),
        -std::log(x)+.5L+x+2*std::log1p(-x),1e-10L,"independent whole-plus convolution");
  close(p.integral(0,.2L)+p.integral(.2L,.7L)+p.integral(.7L,1),0,1e-10L,"plus-bin partition");
  const Distribution delta{1};
  close(delta.integral(.9L,1),1,1e-16L,"endpoint bin includes delta");
  close(delta.integral(0,.9L),0,1e-16L,"interior bin excludes delta");
  close(delta.convolve([](Real z){return z*(1-z);},x),x*(1-x),1e-16L,"delta convolution identity");
  throws<std::invalid_argument>([&]{p.integral(.5L,.5L);},"empty bin");
  throws<std::invalid_argument>([&]{p.convolve({},x);},"empty convolution function");
}
