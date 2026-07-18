// Fixed-Lorentz-component gate (gamma^0 and friends — the finite-T 3+1 split in which the
// four-quark Fierz bases are written), through the FULL NumTracer pipeline.
// Kernels from tests/gen/gen_fixedcomp_numeric.wls.
//
// A fixed component is a Lorentz index pinned to a CONCRETE value. NumTrace rewrites it into a
// contraction with the constant unit basis vector e_i, so gamma^i is emitted as an ordinary slash
// and the C++ engine is untouched. These oracles check the rewrite is numerically faithful.
//
// Frame (propFrameFT, temporal components NONZERO — at a T=0 vacuum frame every q_0 is 0 and every
// oracle below degenerates to 0, i.e. the gate would pass vacuously):
//     q1 = {p0, p, 0, 0}          ql = {l0, l1*cos1, l1*sin1, 0}
#include "Fixc_4_kernel.hh"
#include "Fixc_4ref_kernel.hh"
#include "Fixc_g_kernel.hh"
#include "Fixc_t_kernel.hh"
#include "Fixc_tt_kernel.hh"
#include "Fixc_x_kernel.hh"

#include <cmath>
#include <cstdio>

struct Reg {}; // no regulator/dressing is referenced by these kernels

int main()
{
  struct Pt { double p0, p, l0, l1, cos1; };
  static const Pt pts[] = {
      {0.70, 1.30, -0.40, 0.90, 0.25},  {1.90, 0.55, 1.10, 1.70, -0.60},
      {-0.35, 2.10, 0.80, 0.45, 0.85},  {1.25, 0.95, -1.60, 2.30, -0.15},
      {2.40, 1.75, 0.30, 1.05, 0.70},   {0.15, 0.60, 2.05, 0.80, -0.90},
  };

  double worst = 0.0;
  auto check = [&](const char *what, double got, double want) {
    const double d = std::fabs(got - want);
    worst = std::max(worst, d);
    if (d > 1e-9) std::printf("  MISMATCH %-28s got=%+.12f want=%+.12f\n", what, got, want);
  };

  double worstRef = 0.0;
  for (const Pt &t : pts) {
    const double a[] = {t.p0, t.p, t.l0, t.l1, t.cos1};

    // tr(gamma^0 q1_slash) = 4 * q1_0 = 4 p0   -- and the SPATIAL one, 4 * q1_1 = 4 p.
    // p0 and p are independent args, so an off-by-one in the unit vector swaps these two.
    check("tr(g^0 q1)", numtracer_kernels::Fixc_t_kernel<Reg>::kernel(a[0], a[1], a[2], a[3], a[4]),
          4.0 * t.p0);
    check("tr(g^1 q1)", numtracer_kernels::Fixc_x_kernel<Reg>::kernel(a[0], a[1], a[2], a[3], a[4]),
          4.0 * t.p);

    // tr(gamma^0 gamma^0) = 4 (e_0 . e_0) = 4, a pure constant: the two fixed components must get
    // distinct dummies (else they self-contract) but share one env momentum.
    check("tr(g^0 g^0)", numtracer_kernels::Fixc_tt_kernel<Reg>::kernel(a[0], a[1], a[2], a[3], a[4]),
          4.0);

    // tr(gamma^0 q1_slash gamma^0 ql_slash) = 4[2 q1_0 ql_0 - (q1.ql)] = 4[p0 l0 - p l1 cos1].
    const double q1ql = t.p0 * t.l0 + t.p * t.l1 * t.cos1;
    const double want4 = 4.0 * (2.0 * t.p0 * t.l0 - q1ql);
    const double got4 = numtracer_kernels::Fixc_4_kernel<Reg>::kernel(a[0], a[1], a[2], a[3], a[4]);
    check("tr(g^0 q1 g^0 ql)", got4, want4);

    // The convention-independent cross-check: the SAME trace with gamma^0 written the long way, as
    // gamma^mu contracted against an explicit unit-vector momentum in the frame. That route uses
    // none of the new code path, so agreement proves the rewrite faithful regardless of whether the
    // Euclidean sign conventions in the oracle above are right.
    const double gotRef =
        numtracer_kernels::Fixc_4ref_kernel<Reg>::kernel(a[0], a[1], a[2], a[3], a[4]);
    worstRef = std::max(worstRef, std::fabs(got4 - gotRef));

    // g^{0 nu} q1_nu = q1_0: a fixed component on a METRIC rather than a gamma.
    check("g^{0v} q1_v tr(g^0 ql)",
          numtracer_kernels::Fixc_g_kernel<Reg>::kernel(a[0], a[1], a[2], a[3], a[4]),
          4.0 * t.l0 * t.p0);
  }

  std::printf("fixed-component traces: worst |got - analytic| = %.3e\n", worst);
  std::printf("fixed vs explicit-unit-vector formulation: worst |diff| = %.3e\n", worstRef);

  const bool ok = worst < 1e-9 && worstRef < 1e-9;
  if (!ok)
    std::printf("  (a mismatch on exactly ONE of tr(g^0 q1)/tr(g^1 q1) => the unit basis vector is\n"
                "   off by one component; disagreement with the reference formulation => the\n"
                "   fixed-component rewrite is not equivalent to an explicit e_i contraction)\n");
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
