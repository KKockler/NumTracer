// step-07 — The general 3-point frame reduces to the symmetric-point one.
//
// `frames.wls` generated the SAME network on two different frames:
//
//   sp3_kernel (l1, cos1, cos2, p)                 -- sp3Frame,  symmetric point
//   gen3_kernel(l1, cos1, cos2, p1m, p2m, cosP)    -- gen3Frame, general configuration
//
// The symmetric point is the special case of the general configuration in which the two
// independent external magnitudes are equal and their opening angle is 120 degrees:
//
//   p1m == p2m == p    and    cosP == -1/2 .
//
// This program checks that identity over a grid. It matters because gen3Frame exists precisely to
// escape the symmetric point — some tensor bases are DEGENERATE there (their Gram determinant
// vanishes, so the projectors are not defined). Confirming that the general frame reproduces the
// symmetric one where both are valid is what licenses trusting it where only one is.
#include "gen3_kernel.hh"
#include "sp3_kernel.hh"

#include <cmath>
#include <cstdio>

int main() {
  double worst = 0.0, scale = 0.0;
  int bad = 0;

  // The symmetric point, expressed in gen3Frame's coordinates.
  const double cosP_sym = -0.5;

  for (int ip = 1; ip <= 4; ++ip) {
    const double p = 0.5 * ip;
    for (int il = 1; il <= 4; ++il) {
      const double l1 = 0.5 * il;
      for (int i1 = 0; i1 <= 8; ++i1) {
        const double cos1 = -0.9 + 0.225 * i1;
        for (int i2 = 0; i2 <= 8; ++i2) {
          const double cos2 = -0.9 + 0.225 * i2;

          const double sp = sp3_kernel::kernel(l1, cos1, cos2, p);
          const double gen = gen3_kernel::kernel(l1, cos1, cos2, /*p1m*/ p, /*p2m*/ p, cosP_sym);

          worst = std::max(worst, std::fabs(sp - gen));
          scale = std::max(scale, std::fabs(sp));
          if (std::fabs(sp - gen) >= 1e-11 * std::max(1.0, std::fabs(sp))) ++bad;
        }
      }
    }
  }

  std::printf("sp3Frame  vs  gen3Frame at (p1m = p2m = p, cosP = -1/2)\n");
  std::printf("  1296 kinematic points, worst |difference| = %.3e   (typical |value| ~ %.3g)\n",
              worst, scale);
  std::printf("  points disagreeing beyond 1e-11 relative : %d\n", bad);

  // Off the symmetric point the two must DISAGREE — otherwise the general frame is not general and
  // the test above would be vacuous. This is the guard against a frame that silently ignores its
  // extra arguments.
  const double off = gen3_kernel::kernel(1.0, 0.3, 0.2, /*p1m*/ 1.0, /*p2m*/ 1.7, /*cosP*/ -0.2);
  const double sym = sp3_kernel::kernel(1.0, 0.3, 0.2, 1.0);
  const bool genuinely_general = std::fabs(off - sym) > 1e-6;
  std::printf("  off the symmetric point they differ      : %s (%.6g vs %.6g)\n",
              genuinely_general ? "yes" : "NO — the test above is vacuous!", off, sym);

  const bool ok = (bad == 0) && genuinely_general;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
