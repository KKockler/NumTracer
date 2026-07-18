// Codegen gate for a MULTI-TERM (Plus) Dirac projector of MIXED gamma parity: validate the NUMERIC
// backend kernel generated for the hand-built net in tests/gen/gen_multiterm_proj_numeric.wls
// against a DIRECT ANALYTIC oracle. The net is hand-built (not a FunKit/FORM physics flow), so the
// independent reference is the closed form evaluated from the sp3Frame kinematics by hand — fully
// independent of the matrix-product engine.
//
// The .wls driver already asserts the diagram survives the odd-trace filter at all (it did not
// before the branch-aware parity fix: the diagram-global gamma count was 3 -> odd -> the whole
// diagram was discarded, taking its non-vanishing EVEN branch with it). This file asserts the
// stronger property: the surviving kernel is numerically RIGHT, not merely non-empty.
//
//   sp3Frame(p, l1, cos1, cos2):  q1 = p{1,0,0,0},  q3 = p{-1/2, -sqrt3/2, 0, 0},
//                                 ql = l1{cos1, s1*cos2, s1*sqrt(1-cos2^2), 0},  s1 = sqrt(1-cos1^2)
//
//   proj = 1 + q3_slash                                    (branch gamma-counts 0 and 1)
//   d1   = tr[ proj . q1_slash . ql_slash ]
//        = tr[q1_slash ql_slash] + tr[q3_slash q1_slash ql_slash]
//        = 4 (q1.ql)             + 0                       (the odd-gamma branch vanishes)
//
// so the kernel must equal 4 (q1.ql) exactly. Note this also pins the ODD branch to zero: if the
// 3-gamma term leaked in with any nonzero weight, the q3-dependence would show up as a mismatch.
// (Euclidean metric: a.b = sum_i a_i b_i)
#include "Multiterm_proj_num_kernel.hh"
// Second kernel: a COLLAPSED spinor loop (a closed ntDeltaDirac loop, which orderDiracFacs strips of
// every token) times a colour δ-loop times a Lorentz scalar. Guards compileDirac's nEmptyLoops
// restoration of tr(1) = 4 — the runtime's split_loops discards empty segments, so without it the
// kernel is exactly 4x too small per collapsed loop (the bare-path analogue of a compensation the
// dressed path already had). Oracle: 4 (tr[1]) * 3 (colour Nc) * (q1.ql).
#include "Deltaloop_num_kernel.hh"

#include <array>
#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by this kernel

namespace {
using V4 = std::array<double, 4>;
double dot(const V4 &a, const V4 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }

// analytic oracle for the hand-built multi-term-projector net, over the sp3Frame kinematics
double oracle(double l1, double cos1, double cos2, double p)
{
  const double s1 = std::sqrt(1.0 - cos1 * cos1);
  const V4 q1 = {p, 0.0, 0.0, 0.0};
  const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
  return 4.0 * dot(q1, ql); // tr[q1_slash ql_slash]; the q3_slash branch is an odd trace -> 0
}
} // namespace

int main()
{
  using Num = numtracer_kernels::Multiterm_proj_num_kernel<Reg>;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> Ul(0.05, 3.0), Uc(-0.999, 0.999);
  double pw = 0, maxAbs = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng), cos1 = Uc(rng), cos2 = Uc(rng), p = Ul(rng);
    const double a = oracle(l1, cos1, cos2, p);
    const double n = Num::kernel(l1, cos1, cos2, p);
    pw = std::max(pw, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbs = std::max(maxAbs, std::fabs(a));
  }
  std::printf("multi-term projector codegen numeric vs analytic:  pointwise=%.3e  maxAbs=%.3e\n", pw, maxAbs);

  // collapsed delta-loop: tr(1)=4 must survive a token-free spinor loop
  using Del = numtracer_kernels::Deltaloop_num_kernel<Reg>;
  std::mt19937_64 rng2(90210);
  double pwD = 0, maxAbsD = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng2), cos1 = Uc(rng2), cos2 = Uc(rng2), p = Ul(rng2);
    const double s1 = std::sqrt(1.0 - cos1 * cos1);
    const V4 q1 = {p, 0.0, 0.0, 0.0};
    const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
    const double a = 4.0 * 3.0 * dot(q1, ql); // tr[1] * colour Nc * (q1.ql)
    const double n = Del::kernel(l1, cos1, cos2, p);
    pwD = std::max(pwD, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbsD = std::max(maxAbsD, std::fabs(a));
  }
  std::printf("collapsed delta-loop tr(1)=4 vs analytic:          pointwise=%.3e  maxAbs=%.3e\n", pwD, maxAbsD);
  if (!(pwD < 1e-10 && maxAbsD > 1e-6)) {
    std::printf("TESTS FAILED (collapsed delta-loop: a 0.25 ratio means tr(1)=4 was dropped)\n");
    return 1;
  }
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch — which is exactly
  // the failure mode this test exists for.
  const bool ok = pw < 1e-10 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
