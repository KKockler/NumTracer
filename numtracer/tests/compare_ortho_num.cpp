// Projector-orthonormality gate, through the FULL NumTracer pipeline (kernel from
// tests/gen/gen_ortho_numeric.wls).
//
// A tensor basis's dual projectors must satisfy  <P_i, V_j> = KroneckerDelta[i,j]  exactly and
// INDEPENDENTLY of the kinematics (V_j = TBGetVertex; the dual of a vertex basis is dual to the
// vertex, not to the bare basis element). Rather than commit n^2 kernels, the generator contracts the
// whole matrix against a fixed weight matrix of distinct primes in ONE net:
//
//     W = Sum_ij w[i,j] <P_i, V_j>        must equal        Sum_i w[i,i]
//
// so a single scalar covers every one of the n^2 pairs. The weights are distinct primes, so no
// plausible error cancels: any one wrong entry shifts W by its own weight, and the diagonal weights
// are disjoint from the off-diagonal ones.
//
// AqbqDirect147: n = 3, w[i,j] = Prime[3*(i-1)+j], so the diagonal weights are
// Prime[1] + Prime[5] + Prime[9] = 2 + 11 + 23 = 36.
//
// Traced in a GENERAL 3-point frame (gen3Frame), never the symmetric point: sp3Frame pins a
// measure-zero slice on which the full AqbqDirect basis is degenerate (Det[Gram] = 0, so every dual
// projector is 0/0). Checking orthonormality only there would exercise the tracer exactly where the
// algebra breaks down. The kernels therefore take the general kinematics (p1m, p2m, cosP).
//
// This is the guard that would have caught the multi-term-projector bug on REAL physics bases: that
// bug made <P_7,T_7> trace to an identically-zero kernel, which here would drop W by its weight.
// W must also be momentum-INDEPENDENT — an orthonormality contraction cannot depend on kinematics —
// so any drift across the sampled points is itself a failure.
#include "Ortho147_kernel.hh"
// The FULL AqbqDirect basis: 12 structures, oblique Gram (16 nonzero off-diagonals), all 144 pairs.
// Diagonal weights Prime[12*(i-1)+i], i=1..12, sum to 4543.
#include "Ortho12_kernel.hh"

// SAMPLING. Not random: the Gram of the full basis vanishes on a SURFACE through the symmetric
// point, and uniform sampling lands near it in ~0.03% of draws. There the inverse metric is
// catastrophically ill-conditioned (measured: Det = 39.8 at such a point vs 1.3e20 at a typical one,
// 19 orders down) and the weighted trace loses ~0.3 absolute. That is a conditioning property of the
// BASIS, not a tracer error, but it makes a random-sampling test flaky. So evaluate at a fixed grid
// of well-conditioned points instead — deterministic, and still covering all n^2 pairs per point.
#include <cmath>
#include <cstdio>
#include <random>

struct Reg {}; // no regulator/dressing is referenced by this kernel

int main()
{
  using K = numtracer_kernels::Ortho147_kernel<Reg>;
  const double want = 36.0; // Prime[1] + Prime[5] + Prime[9]

  struct Pt { double l1, cos1, cos2, p1m, p2m, cosP; };
  static const Pt pts[] = {
      {0.70, 0.30, -0.20, 1.00, 1.75, -0.30}, {1.30, -0.55, 0.40, 2.20, 0.60, 0.10},
      {0.40, 0.80, 0.10, 0.50, 2.40, -0.70},  {2.10, -0.20, -0.60, 1.60, 1.10, 0.45},
      {0.90, 0.05, 0.75, 2.80, 2.00, -0.15},  {1.70, -0.85, 0.25, 0.80, 0.35, 0.60},
      {0.25, 0.60, -0.45, 1.20, 2.60, 0.20},  {2.50, 0.15, 0.90, 2.00, 1.40, -0.80},
  };
  double worst = 0.0, spread = 0.0, first = 0.0;
  for (std::size_t i = 0; i < sizeof(pts) / sizeof(pts[0]); ++i) {
    const Pt &t = pts[i];
    const double v = K::kernel(t.l1, t.cos1, t.cos2, t.p1m, t.p2m, t.cosP);
    if (i == 0) first = v;
    worst = std::max(worst, std::fabs(v - want));
    spread = std::max(spread, std::fabs(v - first)); // must be 0: no kinematic dependence
  }
  std::printf("AqbqDirect147 weighted orthonormality trace: got=%.12f want=%.1f\n", first, want);
  std::printf("  worst |W - sum_i w_ii| = %.3e   kinematic spread = %.3e\n", worst, spread);

  using K12 = numtracer_kernels::Ortho12_kernel<Reg>;
  const double want12 = 4543.0; // sum_i Prime[12*(i-1)+i]
  double worst12 = 0.0, spread12 = 0.0, first12 = 0.0;
  for (std::size_t i = 0; i < sizeof(pts) / sizeof(pts[0]); ++i) {
    const Pt &t = pts[i];
    const double v = K12::kernel(t.l1, t.cos1, t.cos2, t.p1m, t.p2m, t.cosP);
    if (i == 0) first12 = v;
    worst12 = std::max(worst12, std::fabs(v - want12));
    spread12 = std::max(spread12, std::fabs(v - first12));
  }
  std::printf("AqbqDirect (full, n=12) weighted trace:      got=%.9f want=%.1f\n", first12, want12);
  std::printf("  worst |W - sum_i w_ii| = %.3e   kinematic spread = %.3e\n", worst12, spread12);

  const bool ok = worst < 1e-9 && spread < 1e-9 &&
                  worst12 < 1e-9 * want12 && spread12 < 1e-9 * want12;
  if (!ok)
    std::printf("  (a deficit near one weight => that pair traced wrong; a nonzero spread => the\n"
                "   contraction picked up kinematics, so it is not an orthonormality relation)\n");
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
