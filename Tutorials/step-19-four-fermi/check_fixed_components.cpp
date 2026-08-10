// step-19 — Fixed-component gammas, checked two independent ways.
//
// `fixed_components.wls` emitted six kernels on the finite-T frame
//
//     q1 = (p0, p, 0, 0),        ql = (l0, l1 c, l1 s, 0),   c = cos1, s = sqrt(1 - c^2)
//
// so q1.ql = p0 l0 + p l1 c, and the temporal components are INDEPENDENT of the spatial ones —
// which is what makes the component-indexing checks below non-vacuous.
//
// Two kinds of check:
//   (1) against analytic oracles, which depend on my Euclidean sign conventions being right;
//   (2) Fx_4 vs Fx_4ref — the SAME trace with gamma^0 written as a fixed component and as an
//       explicit contraction with a unit-vector momentum. That comparison uses none of the new
//       code path on one side, so it is convention-independent: it cannot be fooled by a sign
//       error I also made in (1).
#include "Fx_4_kernel.hh"
#include "Fx_4ref_kernel.hh"
#include "Fx_g_kernel.hh"
#include "Fx_t_kernel.hh"
#include "Fx_tt_kernel.hh"
#include "Fx_x_kernel.hh"

#include <cmath>
#include <cstdio>

namespace {
int failures = 0;
void check(const char *what, double got, double want, const char *form) {
  const double err = std::fabs(got - want) / std::max(1.0, std::fabs(want));
  const bool ok = err < 1e-12;
  if (!ok) ++failures;
  std::printf("  %-12s = %12.6g   expected %12.6g   %-26s %s\n", what, got, want, form,
              ok ? "ok" : "<-- MISMATCH");
}
} // namespace

int main() {
  // Temporal and spatial scales deliberately far apart: an off-by-one in the unit vector would
  // swap p0 and p, and with p0 ~ p the swap could hide inside the tolerance.
  const double p0 = 0.31, p = 2.7, l0 = 0.83, l1 = 1.9, cos1 = 0.4;
  const double q1ql = p0 * l0 + p * l1 * cos1; // q1.ql in this frame

  std::printf("fixed-component gammas at (p0, p, l0, l1, cos1) = (%g, %g, %g, %g, %g)\n", p0, p, l0,
              l1, cos1);

  // B — 0-based component indexing. If e_i landed at slot i+1 these two would swap.
  check("tr(g^0 q1/)", Fx_t_kernel::kernel(p0, p, l0, l1, cos1), 4.0 * p0, "4 q1_0 = 4 p0");
  check("tr(g^1 q1/)", Fx_x_kernel::kernel(p0, p, l0, l1, cos1), 4.0 * p, "4 q1_1 = 4 p");

  // C — a fixed component against itself: 4 (e_0 . e_0) = 4, a pure constant.
  check("tr(g^0 g^0)", Fx_tt_kernel::kernel(p0, p, l0, l1, cos1), 4.0, "4 (e_0.e_0) = 4");

  // D — fixed components interleaved with ordinary slashes.
  check("tr(g^0 q1/ g^0 ql/)", Fx_4_kernel::kernel(p0, p, l0, l1, cos1),
        4.0 * (2.0 * p0 * l0 - q1ql), "4[2 q1_0 ql_0 - q1.ql]");

  // F — a fixed component on a metric rather than a gamma.
  check("g^{0n} q1_n ...", Fx_g_kernel::kernel(p0, p, l0, l1, cos1), 4.0 * l0 * p0, "4 ql_0 p0");

  // E — THE convention-independent check. Fx_4ref writes gamma^0 the long way, as gamma^mu against
  // an explicit unit-vector momentum in the user's frame, and so exercises none of the rewrite.
  std::printf("the convention-independent cross-check\n");
  const double d = Fx_4_kernel::kernel(p0, p, l0, l1, cos1);
  const double e = Fx_4ref_kernel::kernel(p0, p, l0, l1, cos1);
  const bool agree = std::fabs(d - e) <= 1e-12 * std::max(1.0, std::fabs(d));
  if (!agree) ++failures;
  std::printf("  fixed-component  = %12.6g\n  explicit e_0     = %12.6g   %s\n", d, e,
              agree ? "identical" : "<-- MISMATCH");

  std::printf(failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return failures == 0 ? 0 : 1;
}
