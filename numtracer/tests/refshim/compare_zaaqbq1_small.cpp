// NumTracer — the SUB-TERM DEDUP gate (the dense flow).
//
// This is the only flow in the suite with real sub-term density: ~7,500 trace contractions over only
// ~820 distinct traces (9x redundancy), 6 diagrams, and a longest per-net fold of ~140 terms. Every
// other committed flow has exactly ONE sub-term per net, which makes the global sub-term dedup a
// no-op there — the merge, the trace cache and the balanced tree fold are all bypassed, so a green
// suite says nothing about them. Here they all run.
//
// WHAT IS GRADED. The same flow is generated twice by tests/gen/gen_zaaqbq1_small_numeric.wls:
//
//   ZAAqbq1Small_num_kernel     dedup ON  — each distinct trace contracted once, nets fold the table
//   ZAAqbq1SmallRef_num_kernel  dedup OFF — every (net, sub-term) contracted separately (NT_GEN_NO_DEDUP)
//
// The dedup must not change the answer, so the two must agree. Grading is POINTWISE and tight, which
// is deliberate and is NOT the mistake the other compare_* harnesses warn about: those grade against
// a FORM oracle derived independently, where FunKit's FRoute may pick a different loop momentum and
// the integrands then legitimately differ pointwise by O(1) while integrating to the same value (see
// FUNKIT_ROUTING_ISSUE.md). Here BOTH kernels come from the same derivation, the same routing and the
// same Dirac engine — only the codegen scheduling differs — so they must agree term for term, up to
// the reassociation of the tree fold (~1e-15, worst measured 6.7e-15).
//
// A byte-diff of the two headers would NOT work: the dedup changes the order in which traces are
// lowered, hence GlobalEnv's first-seen interning, so every sN renumbers even where the values are
// identical. Judge on numbers, never on bytes.
//
// SCOPE, honestly stated: this validates the DEDUP MACHINERY against the pre-dedup computation. It is
// not an independent physics oracle for the AAqbq vertex — the Dirac engine itself is validated by
// the eleven FORM-oracle'd flows. If the engine mis-traced {[g_mu,g_nu],g_rho}, both sides here would
// be wrong together.
#include "ZAAqbq1Small_num_kernel.hh"    // dedup ON  (gen/)
#include "ZAAqbq1SmallRef_num_kernel.hh" // dedup OFF (gen/) — the control
#include "shim.hpp"

// HOW THE ERROR IS MEASURED, and why not pointwise-relative. This kernel crosses zero: over a random
// sample its magnitude spans ~1e-4 to ~3e+4. At a point that happens to sit near a crossing, ANY last-
// ulp difference is a huge RELATIVE difference — measured, the pointwise relative deviation peaks at
// 3e-11 purely from reassociation, while its MEDIAN is 2.6e-16 (one machine epsilon). So a pointwise
// relative tolerance grades the zero crossings, not the dedup. Two measures that do mean something:
//
//   (1) scale-relative:  max|a-b| / max|kernel|   — the deviation against the kernel's own scale,
//       rather than against whatever tiny value one unlucky sample took. Measured: 8.4e-16, i.e. under
//       four machine epsilons.
//   (2) median pointwise relative deviation — the typical point, where no cancellation is in play.
//       Measured: 2.6e-16.
//
// Both are tight where it counts: a real dedup bug (a dropped trace, a double-counted one, a scalar
// merged onto the wrong key) moves a COEFFICIENT of the polynomial, so it perturbs the kernel by O(1)
// almost everywhere — it would blow (1) and (2) by ten-plus orders of magnitude, not by a factor of a
// few. Loosening a pointwise tolerance to 1e-10 would have hidden nothing here, but it would have been
// measuring the wrong thing; these two measure the right thing and stay at machine precision.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <vector>

int main()
{
  using Dedup = DiFfRG::ZAAqbq1Small_num_kernel;
  using Ref = DiFfRG::ZAAqbq1SmallRef_num_kernel;
  DressingSet d;

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999), Up(0.0, 6.2831853071795865);

  const int Nt = 20000;
  std::vector<double> rel;
  rel.reserve(Nt);
  double maxAbs = 0.0, maxMag = 0.0;

  for (int i = 0; i < Nt; ++i) {
    const double l1 = U(rng), c1 = Uc(rng), c2 = Uc(rng), ph = Up(rng), p = U(rng), k = U(rng);
    // std::real, not a bare double: whether this flow's kernel comes out real- or complex-TYPED
    // depends on which traces survive (the i-bookkeeping — see PruneRealTraces), and that shifted
    // when the generator gained its FMakeSymmetryList. The VALUE is unaffected: measured max |Im|
    // exactly 0 over 20,000 random points, max |Re| 3.5e4. std::real is also a no-op on a double,
    // so this compiles either way — same convention as compare_za4_147_num.cpp.
    const double a = std::real(Dedup::kernel(l1, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAAqbq1,
                                             d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
    const double b = std::real(Ref::kernel(l1, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAAqbq1,
                                           d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));

    if (!std::isfinite(a) || !std::isfinite(b)) {
      std::printf("ZAAqbq1(small): NON-FINITE at l1=%g c1=%g c2=%g phi=%g p=%g k=%g (%g vs %g)\n", l1, c1, c2,
                  ph, p, k, a, b);
      return 1;
    }
    rel.push_back(std::fabs(a - b) / (1e-300 + std::fabs(b)));
    maxAbs = std::max(maxAbs, std::fabs(a - b));
    maxMag = std::max(maxMag, std::fabs(b));
  }

  std::sort(rel.begin(), rel.end());
  const double medRel = rel[rel.size() / 2];
  const double scaleRel = maxAbs / (1e-300 + maxMag);

  std::printf("ZAAqbq1(small) sub-term dedup vs dedup-off control, %d points:\n", Nt);
  std::printf("  scale-relative  max|a-b|/max|k| = %.3e   (kernel scale %.2e)\n", scaleRel, maxMag);
  std::printf("  median pointwise relative dev   = %.3e\n", medRel);
  std::printf("  (max pointwise relative dev     = %.3e — at a zero crossing; not graded)\n", rel.back());

  // ~450x machine epsilon: far above the reassociation floor (both measures sit at 1-4 eps), far below
  // anything a dedup bug can produce.
  constexpr double tol = 1e-13;
  if (!(scaleRel < tol) || !(medRel < tol)) {
    std::printf("  FAIL: dedup changed the kernel (scale-rel %.3e, median-rel %.3e, tol %.1e)\n", scaleRel,
                medRel, tol);
    return 1;
  }
  std::printf("  PASS: dedup is value-preserving to machine precision\n");
  return 0;
}
