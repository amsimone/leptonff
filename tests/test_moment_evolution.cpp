#include "leptonff/moment_evolution.hpp"
#include "check.hpp"

namespace {
void matrix_close(const leptonff::MomentMatrix& a, const leptonff::MomentMatrix& b,
                  leptonff::Real tolerance, const std::string& name) {
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j) close(a[i][j],b[i][j],tolerance,name);
}
}
int main() {
  using namespace leptonff;
  using Mode = QEDCouplingMode;
  MomentEvolutionParameters p;
  for (Mode mode : {Mode::fixed,Mode::one_loop,Mode::two_loop}) {
    p.coupling = mode;
    const auto r = evolve_qed_moment(2,p);
    for (const auto& matrix : {r.evolution,r.fragmentation})
      for (const auto& row : matrix)
        close(row[0]+row[1]+row[2],1,2e-10L,"coupled per-parent momentum");
    if (r.step_refinement_shift > 1e-11L)
      throw std::runtime_error("Moment solver refinement failed");
    const Real a = r.a_initial, L = r.log_scale_ratio, b = 2.0L/3;
    if (mode == Mode::fixed) {
      close(r.a_final,a,1e-15L,"fixed coupling");
      close(r.integral_a,a*L,1e-14L,"fixed coupling integral");
    } else if (mode == Mode::one_loop) {
      close(r.a_final,a/(1-b*a*L),1e-14L,"analytic one-loop coupling");
      close(r.integral_a,-std::log1p(-b*a*L)/b,1e-13L,"analytic one-loop integral");
    } else {
      const Real change = r.a_final-a;
      const Real beta_log = std::log1p(change/(b+a));
      close(change/(b*r.a_final*a)+(beta_log-std::log1p(change/a))/(b*b),L,
            2e-10L,"implicit analytic two-loop coupling");
      close(r.integral_a,(std::log1p(change/a)-beta_log)/b,1e-13L,
            "two-loop integral a");
      close(r.integral_a2,beta_log,1e-13L,"two-loop integral a squared");
    }
    const Real v0 = leading_kernel(Species::muon_minus,Species::muon_minus).moment(2);
    const Real v1 = nlo_valence_kernel().moment(2);
    close(r.evolution[0][0]-r.evolution[0][1],
          std::exp(v0*r.integral_a+v1*r.integral_a2),2e-11L,
          "coupled charge-odd sector vs scalar solution");
  }

  // Noncommuting generators: evolve the high-scale interval with its own
  // correctly transported coupling. Boundary is not reapplied to propagators.
  p.alpha_msbar_at_mu0 = .07L;
  p.coupling = Mode::two_loop;
  const auto full = evolve_qed_moment(3,p);
  auto low = p;
  low.mu_frag = 2;
  const auto first = evolve_qed_moment(3,low);
  auto high = p;
  high.mu_0 = low.mu_frag;
  high.alpha_msbar_at_mu0 = 2*std::acos(-1.0L)*first.a_final;
  const auto second = evolve_qed_moment(3,high);
  matrix_close(full.evolution,multiply_moment_matrices(second.evolution,first.evolution),
               2e-10L,"running matrix composition and ordering");
  auto back = p;
  back.mu_0 = p.mu_frag;
  back.mu_frag = p.mu_0;
  back.alpha_msbar_at_mu0 = 2*std::acos(-1.0L)*full.a_final;
  matrix_close(multiply_moment_matrices(evolve_qed_moment(3,back).evolution,full.evolution),
               identity_moment_matrix(),2e-10L,"backward evolution inverse");

  p = {};
  auto equal = p;
  equal.mu_frag = equal.mu_0;
  const auto identity = evolve_qed_moment(3,equal);
  matrix_close(identity.evolution,identity_moment_matrix(),0,"equal-scale propagator");
  for (unsigned i = 0; i < 3; ++i)
    for (unsigned j = 0; j < 3; ++j)
      close(identity.fragmentation[i][j],(i==j ? 1 : 0)
            +identity.a_initial*identity.boundary1[i][j],1e-15L,"equal-scale boundary");
  p.alpha_msbar_at_mu0 = 0;
  matrix_close(evolve_qed_moment(3,p).fragmentation,identity_moment_matrix(),0,
               "zero coupling identity");

  // A solver re-expansion check, distinct from coefficient formula equality.
  Real previous = 0;
  for (Real a : {.002L,.001L}) {
    p.alpha_msbar_at_mu0 = 2*std::acos(-1.0L)*a;
    const auto r = evolve_qed_moment(3,p);
    Real error = 0;
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned j = 0; j < 3; ++j)
        error = std::max(error,std::abs(r.fragmentation[i][j]-(i==j ? 1 : 0)
            -a*r.coefficient1[i][j]-a*a*r.coefficient2[i][j]));
    if (previous && (previous/error < 7 || previous/error > 9))
      throw std::runtime_error("O(alpha^2) expansion residual is not cubic");
    previous = error;
  }
  throws<std::invalid_argument>([] { evolve_qed_moment(1); },"divergent full number matrix");
  throws<std::invalid_argument>([] {
    MomentEvolutionParameters bad; bad.alpha_msbar_at_mu0 = -1; evolve_qed_moment(2,bad);
  },"invalid MSbar coupling");
  throws<std::invalid_argument>([] {
    MomentEvolutionOptions bad; bad.initial_steps = 0; evolve_qed_moment(2,{},bad);
  },"invalid step controls");
  throws<std::runtime_error>([] {
    MomentEvolutionParameters bad; bad.alpha_msbar_at_mu0 = .9L;
    bad.mu_frag = 1e100L; evolve_qed_moment(2,bad);
  },"unsupported strong running");
  std::cout << "Coupled NLO moments, running, ordering and expansion checks passed\n";
}
