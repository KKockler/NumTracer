// step-08 — Every DSL head, checked against its closed form.
//
// `dsl_heads.wls` generated one micro-kernel per head, all on the same one-angle frame
//
//     q1 = (p, 0, 0, 0),   ql = (l1 c, l1 s, 0, 0),   c = cos1,  s = sqrt(1 - c^2)
//
// so the only invariants are  q1.q1 = p^2,  ql.ql = l1^2,  q1.ql = p l1 c.  Each check below
// writes the expected value in terms of those three and nothing else.
#include "H_dirac2_kernel.hh"
#include "H_dirac4_kernel.hh"
#include "H_diracid_kernel.hh"
#include "H_dress_kernel.hh"
#include "H_longproj_kernel.hh"
#include "H_metric_kernel.hh"
#include "H_sunT_kernel.hh"
#include "H_sunf_kernel.hh"
#include "H_transproj_kernel.hh"

#include <cmath>
#include <cstdio>

namespace {
int failures = 0;

void check(const char *what, double got, double want, const char *form) {
  const double err = std::fabs(got - want) / std::max(1.0, std::fabs(want));
  const bool ok = err < 1e-12;
  if (!ok) ++failures;
  std::printf("  %-14s = %12.6g   expected %12.6g   %-28s %s\n", what, got, want, form,
              ok ? "ok" : "<-- MISMATCH");
}

// A dressing is any callable the consumer supplies. The kernel takes it as `const auto&` and only
// ever calls it — it is opaque to the engine. In a real flow this is a spline interpolator.
struct MyZ {
  double operator()(double x) const { return 1.0 + 0.5 * x * x; }
};
} // namespace

int main() {
  // One representative kinematic point; the closed forms hold for any.
  const double p = 1.7, l1 = 0.9, cos1 = 0.35;
  const double s2 = 1.0 - cos1 * cos1;

  const double q1q1 = p * p;         // q1.q1
  const double qlql = l1 * l1;       // ql.ql
  const double q1ql = p * l1 * cos1; // q1.ql

  std::printf("DSL heads at (p, l1, cos1) = (%g, %g, %g)\n", p, l1, cos1);

  std::printf("Lorentz\n");
  check("ntMetric", H_metric_kernel::kernel(l1, cos1, p), q1q1, "q1.q1");
  check("ntTransProj", H_transproj_kernel::kernel(l1, cos1, p), q1q1 * s2, "q1.q1 - (q1.ql)^2/ql.ql");
  check("ntLongProj", H_longproj_kernel::kernel(l1, cos1, p), q1ql * q1ql / qlql, "(q1.ql)^2 / ql.ql");

  // The two projectors are complementary: P^T + P^L = delta. This is a relation between three
  // INDEPENDENTLY generated kernels, so it checks the engine rather than one closed form.
  const double sum = H_transproj_kernel::kernel(l1, cos1, p) + H_longproj_kernel::kernel(l1, cos1, p);
  check("  P^T + P^L", sum, H_metric_kernel::kernel(l1, cos1, p), "== the metric case");

  std::printf("Dirac\n");
  check("ntDeltaDirac", H_diracid_kernel::kernel(l1, cos1, p), 4.0, "tr(1) = 4");
  check("ntGamma x2", H_dirac2_kernel::kernel(l1, cos1, p), 4.0 * q1ql, "4 (q1.ql)");
  // tr(g^a g^b g^c g^d) A_a B_b C_c D_d = 4[(A.B)(C.D) - (A.C)(B.D) + (A.D)(B.C)]
  // with A = C = q1 and B = D = ql:  4[ (q1.ql)^2 - (q1.q1)(ql.ql) + (q1.ql)^2 ].
  check("ntGamma x4", H_dirac4_kernel::kernel(l1, cos1, p),
        4.0 * (2.0 * q1ql * q1ql - q1q1 * qlql), "4[2(q1.ql)^2 - q1^2 ql^2]");

  std::printf("SU(N)\n");
  // Both SU(N) kernels carry a factor ntSP[q1,q1] = p^2 so that they are not compile-time constants.
  check("ntSUNT", H_sunT_kernel::kernel(l1, cos1, p), 4.0 * q1q1, "tr(T^a T^a) q1^2 = 4 p^2");
  check("ntSUNf", H_sunf_kernel::kernel(l1, cos1, p), 24.0 * q1q1, "f^abc f^abc q1^2 = 24 p^2");

  std::printf("scalars and dressings\n");
  MyZ myZ;
  check("dressing", H_dress_kernel::kernel(l1, cos1, p, myZ), myZ(p) * q1ql * q1q1,
        "myZ(p) (q1.ql) (q1.q1)");

  std::printf(failures == 0 ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return failures == 0 ? 0 : 1;
}
