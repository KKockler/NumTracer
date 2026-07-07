// Codegen gate for the longitudinal/transverse projector leaves: validate the NUMERIC backend
// kernel generated for a hand-built net mixing ntLongProj + ntTransProj (see
// tests/gen/gen_proj_numeric.wls) against a DIRECT ANALYTIC oracle. The net is hand-built (not a
// FunKit/FORM physics flow), so the independent reference is the closed form evaluated from the
// sp3Frame kinematics by hand — fully independent of the matrix-product engine.
//
//   sp3Frame(p, l1, cos1, cos2):  q1 = p{1,0,0,0},  q2 = p{-1/2, sqrt3/2, 0, 0},
//                                 ql = l1{cos1, s1*cos2, s1*sqrt(1-cos2^2), 0},  s1 = sqrt(1-cos1^2)
//   d1 = q1 . P_L(ql) . q1 = (q1.ql)^2 / (ql.ql)                     (longitudinal)
//   d2 = q2 . P_T(ql) . q2 = q2.q2 - (q2.ql)^2 / (ql.ql)            (transverse)
//   kernel value = d1 + d2   (Euclidean metric: a.b = sum_i a_i b_i)
#include "Proj_num_kernel.hh"

#include <array>
#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by these kernels

namespace {
using V4 = std::array<double, 4>;
double dot(const V4 &a, const V4 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }

// analytic oracle for the hand-built projector net, evaluated over the sp3Frame kinematics
double oracle(double l1, double cos1, double cos2, double p)
{
  const double s1 = std::sqrt(1.0 - cos1 * cos1);
  const V4 q1 = {p, 0.0, 0.0, 0.0};
  const V4 q2 = {-p / 2.0, p * std::sqrt(3.0) / 2.0, 0.0, 0.0};
  const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
  const double l2 = dot(ql, ql);
  const double d1 = dot(q1, ql) * dot(q1, ql) / l2;                  // q1 . P_L(ql) . q1
  const double d2 = dot(q2, q2) - dot(q2, ql) * dot(q2, ql) / l2;    // q2 . P_T(ql) . q2
  return d1 + d2;
}
} // namespace

int main()
{
  using Num = numtracer_kernels::Proj_num_kernel<Reg>;

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
  std::printf("projector codegen numeric vs analytic:  pointwise=%.3e  maxAbs=%.3e\n", pw, maxAbs);
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch.
  const bool ok = pw < 1e-10 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
