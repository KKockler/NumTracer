// Codegen gate for the FINITE-T ELECTRIC/MAGNETIC VERTEX SPLIT (see tests/gen/gen_emvertex_numeric.wls).
//
// Three projectors on three DIFFERENT external momenta, contracted through a shared three-index
// core — the structure arXiv:1708.03482 App. B puts under the gluonic vertices, and the one
// qcd-codes/finite_T/YangMills/SP_EM depends on. Two exact identities are checked:
//
//   (a) ORTHOGONALITY   <V[MMM], V[EEM+perm]> = 0        (Pi^M Pi^E = 0, leg by leg)
//   (b) COMPLETENESS    <V[TTT],V[TTT]> - sum_s <V[s],V[s]> = 0   (Pi^M + Pi^E = Pi_perp, plus (a))
//
// Together these are what makes the TensorBases dual projector equal the paper's T_i/<T_i,T_i>, and
// what makes the two-term compact vertex reproduce the pre-split one exactly. If either breaks, the
// dressings extracted downstream are mixtures — and nothing in the physics code would say so.
//
// GRADING. (a) and (b) are ZERO statements, so they are graded ABSOLUTELY against the scale set by
// a third, deliberately non-zero kernel — never pointwise-relative, because two roundings of zero
// disagree by 100% and that has cost time in this repo before.
//
// WHERE THE INDEPENDENCE COMES FROM. There is no dense-oracle twin: MakeNTKernel no longer accepts
// "Backend" -> "Dense". Identity (b) supplies it instead, and more sharply — its left-hand side is
// built from the TRANSVERSE leaf (1/q^2 atoms) and its right-hand side from the electric and
// magnetic leaves (1/|q_vec|^2). Two structurally different lowerings of the same number, required
// to agree exactly.
#include "EMVert_comp_kernel.hh"
#include "EMVert_diag_kernel.hh"
#include "EMVert_off_kernel.hh"

#include <cmath>
#include <cstdio>
#include <random>

int main()
{
  using Diag = numtracer_kernels::EMVert_diag_kernel;
  using Off = numtracer_kernels::EMVert_off_kernel;
  using Comp = numtracer_kernels::EMVert_comp_kernel;

  std::mt19937_64 rng(90210);
  // Matsubara components straddle zero; spatial magnitudes stay away from it so |q_vec|^2 in the
  // magnetic projector never underflows.
  std::uniform_real_distribution<double> Ut(-2.5, 2.5), Up(0.2, 3.0), Uc(-0.99, 0.99),
      Uphi(0.0, 6.283185307179586);

  double maxOff = 0.0, maxComp = 0.0, scale = 0.0;
  for (int i = 0; i < 50000; ++i) {
    const double p01 = Ut(rng), p02 = Ut(rng), p03 = Ut(rng), p = Up(rng), f0 = Ut(rng),
                 l1 = Up(rng), cos1 = Uc(rng), phi = Uphi(rng);

    const double d = Diag::kernel(p01, p02, p03, p, f0, l1, cos1, phi);
    const double o = Off::kernel(p01, p02, p03, p, f0, l1, cos1, phi);
    const double c = Comp::kernel(p01, p02, p03, p, f0, l1, cos1, phi);

    scale = std::fmax(scale, std::fabs(d));
    maxOff = std::fmax(maxOff, std::fabs(o));
    maxComp = std::fmax(maxComp, std::fabs(c));
  }

  std::printf("E/M vertex split:\n");
  std::printf("  scale |<V[MMM],V[MMM]>+...| : %.3e\n", scale);
  std::printf("  (a) orthogonality  max|.|   : %.3e\n", maxOff);
  std::printf("  (b) completeness   max|.|   : %.3e\n", maxComp);

  // The zeros are only meaningful relative to a scale that is genuinely large; a kernel that
  // returned zero everywhere would satisfy (a) and (b) for free.
  const bool nonvacuous = scale > 1e-3;
  const bool ok = nonvacuous && maxOff < 1e-8 * scale && maxComp < 1e-8 * scale;
  if (!nonvacuous)
    std::printf("  !! the control kernel is ~0 everywhere: the zero-identities above prove nothing\n");
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
