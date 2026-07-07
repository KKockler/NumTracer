// Codegen gate for a DISCONNECTED Dirac x Lorentz diagram: validate the NUMERIC backend kernel
// generated for a hand-built net whose diagrams each factorize into a closed Dirac trace, a constant
// colour trace, and an index-disjoint pure-Lorentz scalar (see tests/gen/gen_discdirac_numeric.wls)
// against a DIRECT ANALYTIC oracle. This is the regression guard for the MakeNTKernel::disconnectmix
// case (the hSigL sigma-quark vertex, QCD.nb): the numeric backend FACTORS the disconnected Lorentz
// scalar as its own trace and multiplies it in at assembly (coeff * tr_dirac * tr_colour * tr_lorentz).
// The net is hand-built (no FunKit/FORM), so the independent reference is the closed form evaluated
// from the sp3Frame kinematics by hand:
//   d1 = tr(q1 ql) * Nc * (q3 . P_T(ql) . q3) = 4 (q1.ql) * Nc * (q3.q3 - (q3.ql)^2/ql.ql)
//   d2 = tr(q2 ql) * (q1 . P_L(ql) . q1)      = 4 (q2.ql) * (q1.ql)^2/(ql.ql)
//   d3 = tr(q1 ql) * tr(q2 q3)                = 16 (q1.ql)(q2.q3)     [two disconnected Dirac traces]
//   kernel value = d1 + d2 + d3   (Euclidean, tr(slash a slash b) = 4 a.b, colour Tr(1) = Nc = 3)
#include "Discdirac_num_kernel.hh"

#include <array>
#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by these kernels

namespace {
constexpr double Nc = 3.0;
using V4 = std::array<double, 4>;
double dot(const V4 &a, const V4 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }

double oracle(double l1, double cos1, double cos2, double p)
{
  const double s1 = std::sqrt(1.0 - cos1 * cos1), s3 = std::sqrt(3.0);
  const V4 q1 = {p, 0.0, 0.0, 0.0};
  const V4 q2 = {-p / 2.0, p * s3 / 2.0, 0.0, 0.0};
  const V4 q3 = {-p / 2.0, -p * s3 / 2.0, 0.0, 0.0};
  const V4 ql = {l1 * cos1, l1 * s1 * cos2, l1 * s1 * std::sqrt(1.0 - cos2 * cos2), 0.0};
  const double l2 = dot(ql, ql);
  const double d1 = 4.0 * dot(q1, ql) * Nc * (dot(q3, q3) - dot(q3, ql) * dot(q3, ql) / l2);
  const double d2 = 4.0 * dot(q2, ql) * (dot(q1, ql) * dot(q1, ql) / l2);
  const double d3 = 16.0 * dot(q1, ql) * dot(q2, q3);
  return d1 + d2 + d3;
}
} // namespace

int main()
{
  using Num = numtracer_kernels::Discdirac_num_kernel<Reg>;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> Ul(0.05, 3.0), Uc(-0.999, 0.999);
  // This kernel is a SUM of products of independent traces (d1 + d2 + d3); the summands can nearly
  // cancel, so a pure pointwise-relative metric blows up at zeros of the sum even when numeric and
  // the oracle agree to machine precision. Gate on the max ABSOLUTE error relative to the kernel's
  // overall scale (maxAbs) — the standard, cancellation-robust check.
  double maxAbsErr = 0, maxAbs = 0;
  for (int i = 0; i < 200000; ++i) {
    const double l1 = Ul(rng), cos1 = Uc(rng), cos2 = Uc(rng), p = Ul(rng);
    const double a = oracle(l1, cos1, cos2, p);
    const double n = Num::kernel(l1, cos1, cos2, p);
    maxAbsErr = std::max(maxAbsErr, std::fabs(n - a));
    maxAbs = std::max(maxAbs, std::fabs(a));
  }
  const double rel = maxAbsErr / (1e-300 + maxAbs);
  std::printf("disconnected Dirac x Lorentz codegen numeric vs analytic:  maxAbsErr=%.3e  maxAbs=%.3e  scaled=%.3e\n",
              maxAbsErr, maxAbs, rel);
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch.
  const bool ok = rel < 1e-12 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
