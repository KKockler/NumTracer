// Codegen gate for FINITE-TEMPERATURE support: validate the NUMERIC backend kernel generated for a
// finite-T net (see tests/gen/gen_ftproj_numeric.wls) — electric/magnetic/transverse projectors plus
// the spatial scalar product ntSPS and the temporal component ntVec[q,0], over the finite-T frame
// propFrameFT — against a DIRECT ANALYTIC oracle. The net is hand-built (not a FunKit/FORM physics
// flow), so the independent reference is the closed form evaluated from the frame kinematics by hand,
// with the engine's finite-T projector conventions:
//   P_T (4D transverse):        pp.P_T(ll).pp = pp.pp - (pp.ll)^2/(ll.ll)               [4D dots]
//   P_M (spatial magnetic):     pp.P_M(ll).pp = ppS.ppS - (ppS.llS)^2/|llS|^2           [spatial dots]
//   P_E = P_T - P_M (electric): pp.P_E(ll).pp = pp.P_T(ll).pp - pp.P_M(ll).pp
// net = pp.P_E(ll).pp + pp.P_M(ll).pp + sps(pp,ll)*ll_0*(pp.P_T(ll).pp)
//   propFrameFT(p0,p,l0,l1,cos1):  pp = {p0, p, 0, 0},  ll = {l0, l1*cos1, l1*sqrt(1-cos1^2), 0}
//   sps = spatial dot (components 1..3), ll_0 = temporal component. Euclidean metric.
#include "FTProj_num_kernel.hh"

#include <array>
#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by these kernels

namespace {
using V4 = std::array<double, 4>;
double dot4(const V4 &a, const V4 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }
// spatial dot (components 1..3), component 0 is the temporal/Matsubara direction
double dotS(const V4 &a, const V4 &b) { return a[1] * b[1] + a[2] * b[2] + a[3] * b[3]; }

double oracle(double p0, double p, double l0, double l1, double cos1)
{
  const double s1 = std::sqrt(1.0 - cos1 * cos1);
  const V4 pp = {p0, p, 0.0, 0.0};
  const V4 ll = {l0, l1 * cos1, l1 * s1, 0.0};
  const double pT = dot4(pp, pp) - dot4(pp, ll) * dot4(pp, ll) / dot4(ll, ll);        // pp.P_T.pp
  const double pM = dotS(pp, pp) - dotS(pp, ll) * dotS(pp, ll) / dotS(ll, ll);        // pp.P_M.pp
  const double d1 = pT - pM;                                                          // pp.P_E.pp
  const double d2 = pM;                                                               // pp.P_M.pp
  const double d3 = dotS(pp, ll) * ll[0] * pT;                                         // sps*ll_0*pp.P_T.pp
  return d1 + d2 + d3;
}
} // namespace

int main()
{
  using Num = numtracer_kernels::FTProj_num_kernel<Reg>;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> Up(0.05, 3.0), Uc(-0.999, 0.999), Ut(-3.0, 3.0);
  double pw = 0, maxAbs = 0;
  for (int i = 0; i < 200000; ++i) {
    const double p0 = Ut(rng), p = Up(rng), l0 = Ut(rng), l1 = Up(rng), cos1 = Uc(rng);
    const double a = oracle(p0, p, l0, l1, cos1);
    const double n = Num::kernel(p0, p, l0, l1, cos1);
    pw = std::max(pw, std::fabs(n - a) / (1e-300 + std::fabs(a)));
    maxAbs = std::max(maxAbs, std::fabs(a));
  }
  std::printf("finite-T projector codegen numeric vs analytic:  pointwise=%.3e  maxAbs=%.3e\n", pw, maxAbs);
  // maxAbs guards against a degenerate all-zero kernel masking a real mismatch.
  const bool ok = pw < 1e-10 && maxAbs > 1e-6;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
