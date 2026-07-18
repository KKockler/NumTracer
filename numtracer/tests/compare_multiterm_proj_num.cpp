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
// Third kernel: a momentum-free pure-LORENTZ component (a closed metric loop g_{uv} g^{uv} = D)
// times a Lorentz scalar. Guards the "Constant" flag meaning "a constant SU(N) component": a
// momentum-free metric/epsilon component is not colour, and was emitted as raw Mathematica
// (colFacG[ntMetric[...]]) into the generator. Oracle: 4 (= D) * (q1.ql).
#include "Metricloop_num_kernel.hh"
// A CONSTANT colour/flavour component that is a Plus of SU(N) heads — compileColG had no Plus case.
#include "Colsum_num_kernel.hh"
// SU(N) FUNDAMENTAL Levi-Civita: epsilon PAIRS contracted into Kronecker deltas.
#include "Epsfund_num_kernel.hh"
// ADJOINT Levi-Civita -- SU(2) ONLY (eps^abc = f^abc there, and ONLY there).
#include "Epsadj_num_kernel.hh"

#include "Epsstraddle_num_kernel.hh"

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

  // momentum-free metric loop: g_{uv} g^{uv} = D = 4
  using Met = numtracer_kernels::Metricloop_num_kernel<Reg>;
  std::mt19937_64 rng3(5150);
  double pwM = 0, maxAbsM = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng3), cos1 = Uc(rng3), cos2 = Uc(rng3), p = Ul(rng3);
    const double s1 = std::sqrt(1.0 - cos1 * cos1);
    const V4 q1 = {p, 0.0, 0.0, 0.0};
    const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
    const double a = 4.0 * dot(q1, ql); // D * (q1.ql)
    const double n = Met::kernel(l1, cos1, cos2, p);
    pwM = std::max(pwM, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbsM = std::max(maxAbsM, std::fabs(a));
  }
  std::printf("momentum-free metric loop vs analytic:             pointwise=%.3e  maxAbs=%.3e\n", pwM, maxAbsM);
  if (!(pwM < 1e-10 && maxAbsM > 1e-6)) {
    std::printf("TESTS FAILED (metric loop)\n");
    return 1;
  }
  // CONSTANT COLOUR/FLAVOUR SUM: a Plus of SU(N) heads, which compileColG had no case for — it fell
  // through to colFacG (the one per-head dispatcher with neither a Plus branch nor a catch-all) and
  // emitted raw Mathematica into the generator .cpp. Found on the four-quark Fierz flavour structure
  // delta*delta - 4*T*T; reproduced here without FunKit.
  //   delta_{ij} delta_{ji} = Nc = 3,  T^a_{ij} T^a_{ji} = tr(T^a T^a) = (Nc^2-1)/2 = 4
  //   => 3 + (-2)(4) = -5, times the Lorentz factor (q1.ql).
  // The sign matters: a fix that dropped a branch, or summed the branches as a product, would land
  // on +3, -8 or -24 rather than -5, all of which this catches.
  using Csum = numtracer_kernels::Colsum_num_kernel<Reg>;
  std::mt19937_64 rng4(90210);
  double pwC = 0, maxAbsC = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng4), cos1 = Uc(rng4), cos2 = Uc(rng4), p = Ul(rng4);
    const double s1 = std::sqrt(1.0 - cos1 * cos1);
    const V4 q1 = {p, 0.0, 0.0, 0.0};
    const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
    const double a = -5.0 * dot(q1, ql);
    const double n = Csum::kernel(l1, cos1, cos2, p);
    pwC = std::max(pwC, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbsC = std::max(maxAbsC, std::fabs(a));
  }
  std::printf("constant colour SUM (3 - 2*4 = -5) vs analytic:    pointwise=%.3e  maxAbs=%.3e\n", pwC, maxAbsC);
  if (!(pwC < 1e-10 && maxAbsC > 1e-6)) {
    std::printf("TESTS FAILED (constant colour sum: check the compileColGSum branch list)\n");
    return 1;
  }

  // SU(N) FUNDAMENTAL Levi-Civita. NumTrace contracts epsilon PAIRS into Kronecker deltas
  // (eps.eps = k! det(delta)) — a rewrite into existing primitives, so the C++ engine never sees an
  // N-index object. Four closed contractions, distinct-prime weighted:
  //   H1 SU(3) k=1 direct  = 3*3 - 3 =  6      H2 SU(3) k=1 crossed = 3 - 3*3 = -6
  //   H3 SU(2) k=0 direct  = 2*2 - 2 =  2      H4 SU(2) k=0 crossed = 2 - 2*2 = -2
  //   2*6 + 3*(-6) + 5*2 + 7*(-2) = -10
  // Direct vs crossed differ ONLY by the antisymmetry sign, so a dropped Signature flips H2/H4 and
  // lands on +22 rather than perturbing -10 — a sign error cannot hide here.
  using Eps = numtracer_kernels::Epsfund_num_kernel<Reg>;
  std::mt19937_64 rng5(31337);
  double pwE = 0, maxAbsE = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng5), cos1 = Uc(rng5), cos2 = Uc(rng5), p = Ul(rng5);
    const double s1 = std::sqrt(1.0 - cos1 * cos1);
    const V4 q1 = {p, 0.0, 0.0, 0.0};
    const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
    const double a = -10.0 * dot(q1, ql);
    const double n = Eps::kernel(l1, cos1, cos2, p);
    pwE = std::max(pwE, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbsE = std::max(maxAbsE, std::fabs(a));
  }
  std::printf("fundamental epsilon pairs (-10) vs analytic:       pointwise=%.3e  maxAbs=%.3e\n", pwE, maxAbsE);
  if (!(pwE < 1e-10 && maxAbsE > 1e-6)) {
    std::printf("TESTS FAILED (epsilon: +22 instead of -10 means a dropped antisymmetry sign)\n");
    return 1;
  }

  // ADJOINT Levi-Civita. !! SU(2) ONLY !! At rank 2 -- and only at rank 2 -- eps^abc = f^abc
  // exactly (coefficient +1, verified independently from Pauli matrices with T^a = sigma^a/2 and
  // f = -2i tr([T^a,T^b]T^c), matching sun_net.hpp:225,233-247). For SU(3) f is not an epsilon at
  // all (f^123 = 1 but f^147 = 1/2), and rank != 2 is refused by FromFunKit::epsadj.
  //
  // Oracle: eps^abc f^abc = f^abc f^abc = 3! = 6 for SU(2), times (q1.ql). This contraction is
  // LINEAR in the rewrite coefficient, which is the point: every quadratic-in-eps contraction is
  // invariant under c -> -c and would pass with the sign inverted.
  using Eadj = numtracer_kernels::Epsadj_num_kernel<Reg>;
  std::mt19937_64 rng6(24601);
  double pwA = 0, maxAbsA = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng6), cos1 = Uc(rng6), cos2 = Uc(rng6), p = Ul(rng6);
    const double s1 = std::sqrt(1.0 - cos1 * cos1);
    const V4 q1 = {p, 0.0, 0.0, 0.0};
    const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
    const double a = 6.0 * dot(q1, ql);
    const double n = Eadj::kernel(l1, cos1, cos2, p);
    pwA = std::max(pwA, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbsA = std::max(maxAbsA, std::fabs(a));
  }
  std::printf("SU(2) adjoint epsilon = f (6) vs analytic:         pointwise=%.3e  maxAbs=%.3e\n", pwA, maxAbsA);
  if (!(pwA < 1e-10 && maxAbsA > 1e-6)) {
    std::printf("TESTS FAILED (adjoint eps: -6 means the eps->f rewrite flipped sign)\n");
    return 1;
  }

  // Epsilon pair STRADDLING an eager Plus, ranks interleaved -- the diquark shape of the four-quark
  // Fierz bases. The outer epsilons sit outside the sum and their partners inside it, so each
  // summand alone holds one rank-3 and one rank-2 epsilon. A rank-blind pairing joins those two,
  // sizes the determinant by one of them and silently drops a colour leg; the summands then expose
  // different free indices and the run dies downstream in checkLabels (NumTrace::plusfree).
  //
  // Oracle: summand 1 = 2 * (direct colour 6) * (direct flavour 2) = 24
  //         summand 2 = 3 * (crossed colour -6) * (crossed flavour -2) = 36   -> 60 * (q1.ql)
  // The crossed summand carries two sign flips, so one dropped Signature lands on -12, not near 60.
  using Estr = numtracer_kernels::Epsstraddle_num_kernel<Reg>;
  std::mt19937_64 rng7(90210);
  double pwS = 0, maxAbsS = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng7), cos1 = Uc(rng7), cos2 = Uc(rng7), p = Ul(rng7);
    const double s1 = std::sqrt(1.0 - cos1 * cos1);
    const V4 q1 = {p, 0.0, 0.0, 0.0};
    const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
    const double a = 60.0 * dot(q1, ql);
    const double n = Estr::kernel(l1, cos1, cos2, p);
    pwS = std::max(pwS, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbsS = std::max(maxAbsS, std::fabs(a));
  }
  std::printf("epsilon pair straddling a Plus (60) vs analytic:   pointwise=%.3e  maxAbs=%.3e\n", pwS, maxAbsS);
  if (!(pwS < 1e-10 && maxAbsS > 1e-6)) {
    std::printf("TESTS FAILED (straddle: a cross-rank mispairing drops an index and mis-contracts)\n");
    return 1;
  }

  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch — which is exactly
  // the failure mode this test exists for.
  const bool ok = pw < 1e-10 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
