// End-to-end check of the ASSEMBLED quark-gluon vertex flow (A qb q), projected onto EACH external
// tensor structure {1,4,7} of the AqbqDirect147 basis: NumTracer's numeric backend vs the
// FormTracer oracles.
//
// Graded on the ANGLE-INTEGRATED value, not pointwise. The two backends may route the loop momentum
// differently; a routing difference is an odd term that integrates to zero, so it can show up
// pointwise but MUST cancel in the integral. Only the integrated number is physical. (The pointwise
// figure is still printed, as a diagnostic for whether the routings happen to coincide.)
//
// Both sides must be generated from the SAME qcd147_setup.m vintage -- a mixed-vintage pair
// disagrees even integrated. Regenerate with:
//     wolfram -script tests/gen/gen_qcd_aqbq147_form.wls      (writes refshim/*_147_kernel.hh)
//     rm -f tests/gen/ZAqbq*_147_num_kernel.hh                (MakeNTKernel SKIPS existing wrappers)
//     wolfram -script tests/gen/gen_qcd_aqbq147_numeric.wls
//
// Measure: in 4D, dOmega_3 ∝ sin^2(th1) dth1 sin(th2) dth2 -> sqrt(1-c1^2) dc1 · dc2.
// Grading is SCALE-relative (divided by the largest |I| over the sample), never pointwise-relative,
// so a near-zero integral cannot manufacture a huge ratio.
#include "ZAqbq1_147_kernel.hh"     // FormTracer oracle, struct 1 (refshim/)
#include "ZAqbq4_147_kernel.hh"     // FormTracer oracle, struct 4 (refshim/)
#include "ZAqbq7_147_kernel.hh"     // FormTracer oracle, struct 7 (refshim/)
#include "ZAqbq1_147_num_kernel.hh" // numeric backend, struct 1 (gen/)
#include "ZAqbq4_147_num_kernel.hh" // numeric backend, struct 4 (gen/)
#include "ZAqbq7_147_num_kernel.hh" // numeric backend, struct 7 (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

// Gauss-Legendre nodes/weights on [-1,1] via Newton iteration on P_n.
static void gauss_legendre(int n, std::vector<double> &x, std::vector<double> &w)
{
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; ++i) {
    double z = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 0;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1, p1 = 0;
      for (int j = 0; j < n; ++j) {
        double p2 = p1;
        p1 = p0;
        p0 = ((2.0 * j + 1.0) * z * p1 - j * p2) / (j + 1);
      }
      dp = n * (z * p0 - p1) / (z * z - 1.0);
      double dz = -p0 / dp;
      z += dz;
      if (std::fabs(dz) < 1e-15) break;
    }
    x[i] = z;
    w[i] = 2.0 / ((1.0 - z * z) * dp * dp);
  }
}

int main()
{
  DressingSet d;
  const double Nf = 2.0;

  const int NQ = 48;
  std::vector<double> xc, wc;
  gauss_legendre(NQ, xc, wc);

  auto run = [&](int lab, auto formK, auto numK) {
    auto cForm = [&](double l1, double c1, double c2, double p, double k) {
      return formK(l1, c1, c2, p, k, Nf, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7, d.dtZc, d.Zc, d.dtZA,
                   d.ZA, d.dtZq, d.Zq, d.Mq);
    };
    auto cNum = [&](double l1, double c1, double c2, double p, double k) {
      return std::real(numK(l1, c1, c2, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.ZAqbq4, d.ZAqbq7, d.dtZc, d.Zc,
                            d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
    };

    std::mt19937_64 rng(31337);
    std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
    double pw = 0;
    for (int i = 0; i < 100000; ++i) {
      double l1 = U(rng), c1 = Uc(rng), c2 = Uc(rng), p = U(rng), k = U(rng);
      double r = cForm(l1, c1, c2, p, k);
      pw = std::max(pw, std::fabs(cNum(l1, c1, c2, p, k) - r) / (1e-300 + std::fabs(r)));
    }

    std::mt19937_64 rng2(20260719);
    double maxdiff = 0, scale = 0;
    for (int s = 0; s < 300; ++s) {
      double l1 = U(rng2), p = U(rng2), k = U(rng2);
      double IF = 0, IN = 0;
      for (int i = 0; i < NQ; ++i) {
        const double c1 = xc[i], w1 = wc[i] * std::sqrt(1.0 - c1 * c1);
        for (int j = 0; j < NQ; ++j) {
          const double w = w1 * wc[j];
          IF += w * cForm(l1, c1, xc[j], p, k);
          IN += w * cNum(l1, c1, xc[j], p, k);
        }
      }
      scale = std::max(scale, std::fabs(IF));
      maxdiff = std::max(maxdiff, std::fabs(IN - IF));
    }
    const double rel = maxdiff / (1e-300 + scale);
    const bool ok = rel < 1e-10;
    std::printf("  ZAqbq%d (147):  pointwise=%.3e   2-angle-INTEGRATED scale-rel=%.3e  (|scale|=%.2e)  %s\n", lab, pw,
                rel, scale, ok ? "ok" : "FAIL");
    return ok;
  };

  std::printf("QCD quark-gluon vertex flow, {1,4,7} basis -- assembled flow vs FormTracer oracle:\n");
  bool ok = true;
  ok &= run(1, &DiFfRG::ZAqbq1_147_kernel<DiFfRG::ShimRegulator>::kernel,
            &DiFfRG::ZAqbq1_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  ok &= run(4, &DiFfRG::ZAqbq4_147_kernel<DiFfRG::ShimRegulator>::kernel,
            &DiFfRG::ZAqbq4_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  ok &= run(7, &DiFfRG::ZAqbq7_147_kernel<DiFfRG::ShimRegulator>::kernel,
            &DiFfRG::ZAqbq7_147_num_kernel<DiFfRG::ShimRegulator>::kernel);
  std::printf(ok ? "\nINTEGRATED CHECK PASSED\n" : "\nINTEGRATED CHECK FAILED\n");
  return ok ? 0 : 1;
}
