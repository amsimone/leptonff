#include "leptonff/distribution.hpp"
#include "check.hpp"

#include <limits>

int main() {
  using namespace leptonff;
  const Real pi = std::acos(-1.0L);
  const Distribution delta(2);
  const Distribution regular(0, [](Real z) { return 3*z*z; });
  const Distribution soft(0, {}, {[](Real z) { return 1/(1-z); }});
  const Distribution logsoft(0, {}, {[](Real z) { return std::log1p(-z)/(1-z); }});
  close(delta.moment(3), 2, 1e-15L, "delta moment");
  close(delta.resolved_density(.4L),0,0,"delta has no interior density");
  close((regular+soft).resolved_density(.4L),3*.4L*.4L+1/.6L,1e-14L,
        "interior density excludes endpoint and retains plus kernel");
  throws<std::invalid_argument>([&]{soft.resolved_density(1);},"interior excludes plus endpoint");
  close(regular.moment(2), .75L, 1e-12L, "regular moment");
  close(soft.moment(1), 0, 1e-15L, "plus on a constant");
  close(soft.moment(4), -11.0L/6, 1e-10L, "plus harmonic moment");
  close(logsoft.moment(2), 1, 1e-9L, "log-plus moment");
  close(soft.apply([](Real) { return 1; }, .3L), std::log(.7L), 1e-10L,
        "finite lower bound endpoint");
  close(soft.apply([](Real z) { return z; }, .3L), -.7L + std::log(.7L),
        1e-10L, "finite lower bound nonconstant");
  close(soft.apply([](Real z) { return z > .8L ? 1 : 0; }, 0, {}, {.8L}),
        std::log(.2L), 1e-10L, "step measurement with breakpoint");
  close((delta + regular + soft.scaled(-2)).moment(2), 4.75L, 1e-10L,
        "signed linear combination");

  const auto coupled = [](Real z1, Real z2) { return z1+z2 > 1 ? 1.0L : 0.0L; };
  const auto breaks = [](Real z1) { return std::vector<Real>{1-z1}; };
  close(apply_pair(soft, soft, coupled, {}, {}, breaks), -pi*pi/6, 1e-8L,
        "coupled cut double-plus including endpoints");
  const auto mixed = Distribution(1) + soft.scaled(.2L);
  close(apply_pair(mixed, mixed, coupled, {}, {}, breaks), 1-.04L*pi*pi/6,
        1e-8L, "delta plus cross terms");
  close(apply_pair(soft, logsoft, [](Real x, Real y) { return x*y; }), -1,
        1e-8L, "tensor product moments");
  close(apply_pair_order1(soft, soft, .2L, coupled), 1, 1e-12L,
        "pair order1 excludes quadratic term");

  close(integrate([](Real z) { return std::log(z); }, 0, 1).value, -1,
        1e-9L, "integrable endpoint logarithm");
  close(integrate([](Real z) { return z*z; }, 0, 1, {}, {.7L,.2L,.7L}).value,
        1.0L/3, 1e-12L, "sorted unique breakpoints");
  throws<std::invalid_argument>([&] { soft.apply([](Real) { return 1; }, 1); },
                                "invalid distribution interval");
  throws<std::invalid_argument>([&] { soft.moment(0); }, "invalid moment");
  throws<std::runtime_error>([] {
    integrate([](Real) { return std::numeric_limits<Real>::quiet_NaN(); }, 0, 1);
  }, "nonfinite quadrature");
  throws<std::runtime_error>([] {
    integrate([](Real x) { return std::log(x); }, 0, 1, {1e-15L, 0, 1});
  }, "unconverged quadrature must fail");
  throws<std::invalid_argument>([] {
    integrate([](Real x) { return x; }, 0, 1, {-1, 0, 10});
  }, "invalid tolerance");
  std::cout << "Distribution and quadrature checks passed\n";
}
