// NumTracer — the GENERAL-FRAME (polyFrameSpec) gate.
//
// This is the only flow in the suite off the unit-loop fast path: its externals are parametrised by
// the (S0,S1,SPhi) shape coordinates (the with_mesons lambda1L3D class), so unitLoopOkQ fails and
// the general polyFrameSpec parametrisation runs — minted radical symbols, nsym 12, momentum bases
// in exact +/- pairs, dressed+complex traces. Every mechanism the frame/sign levers touch lives on
// this path, and nothing else covers it. See tests/qcd_lambda3d_setup.m for the physics content.
//
// WHAT IS GRADED. The same flow is generated twice by tests/gen/gen_lambda3d_small_numeric.wls:
//
//   Lambda3DSmall_num_kernel     dedup ON  (tests/gen/lambda3d/)
//   Lambda3DSmallRef_num_kernel  dedup OFF (tests/gen/lambda3d_ref/, NT_GEN_NO_DEDUP) — the control
//
// Both come from the SAME derivation, routing and Dirac engine, so they must agree pointwise up to
// the tree fold's reassociation. As frame/sign levers land, their escape-hatch controls are graded
// through this same harness. A byte-diff would not work — the dedup renumbers every sN; and a
// pointwise-relative tolerance would grade the kernel's zero crossings, not the machinery. The two
// meaningful measures (see compare_zaaqbq1_small.cpp for the full argument):
//   (1) scale-relative max|a-b| / max|kernel|, and (2) the median pointwise relative deviation.
#include "Lambda3DSmall_num_kernel.hh"    // dedup ON  (gen/lambda3d/)
#include "Lambda3DSmallRef_num_kernel.hh" // dedup OFF (gen/lambda3d_ref/) — the control
#include "shim.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

int main()
{
  using Dedup = DiFfRG::Lambda3DSmall_num_kernel;
  using Ref = DiFfRG::Lambda3DSmallRef_num_kernel;
  DressingSet d;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  // S1 stays inside the open unit disc (the triangle/Kaellen condition); SPhi covers the circle.
  std::uniform_real_distribution<double> Us1(0.001, 0.9999), Uphi(-3.1415926535897931, 3.1415926535897931);

  const int Nt = 20000;
  std::vector<double> rel;
  rel.reserve(Nt);
  double maxAbs = 0.0, maxMag = 0.0;

  for (int i = 0; i < Nt; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), c2 = Uc(rng), S0 = U(rng), S1 = Us1(rng), SPhi = Uphi(rng),
                 k = U(rng);
    const double a = Dedup::kernel(l1, c1, c2, S0, S1, SPhi, k, d.ZAqbq1, d.Zq, d.Mq, d.ZA, d.dtZA, d.dtZq);
    const double b = Ref::kernel(l1, c1, c2, S0, S1, SPhi, k, d.ZAqbq1, d.Zq, d.Mq, d.ZA, d.dtZA, d.dtZq);

    if (!std::isfinite(a) || !std::isfinite(b)) {
      std::printf("lambda3d(small): NON-FINITE at l1=%g c1=%g c2=%g S0=%g S1=%g SPhi=%g k=%g (%g vs %g)\n",
                  l1, c1, c2, S0, S1, SPhi, k, a, b);
      return 1;
    }
    rel.push_back(std::fabs(a - b) / (1e-300 + std::fabs(b)));
    maxAbs = std::max(maxAbs, std::fabs(a - b));
    maxMag = std::max(maxMag, std::fabs(b));
  }

  std::sort(rel.begin(), rel.end());
  const double medRel = rel[rel.size() / 2];
  const double scaleRel = maxAbs / (1e-300 + maxMag);

  std::printf("lambda3d(small) polyFrameSpec flow, dedup vs dedup-off control, %d points:\n", Nt);
  std::printf("  scale-relative  max|a-b|/max|k| = %.3e   (kernel scale %.2e)\n", scaleRel, maxMag);
  std::printf("  median pointwise relative dev   = %.3e\n", medRel);
  std::printf("  (max pointwise relative dev     = %.3e — at a zero crossing; not graded)\n", rel.back());

  constexpr double tol = 1e-13;
  if (!(scaleRel < tol) || !(medRel < tol)) {
    std::printf("  FAIL: kernels disagree (scale-rel %.3e, median-rel %.3e, tol %.1e)\n", scaleRel, medRel,
                tol);
    return 1;
  }
  std::printf("  PASS: value-preserving to machine precision\n");
  return 0;
}
