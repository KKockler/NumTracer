// Validate the NUMERIC (matrix-product) backend for the QCD ZA gluon self-energy against the
// independent FormTracer (FORM) oracle — the two tensor-trace engines computing the SAME gluon
// self-energy loop (1 angle). As for the other vacuum flows the pointwise value can differ by a
// loop-routing term odd in cos1 that integrates to zero, so the physical check is the cos1-integrated
// value.
#include "ZA_form_kernel.hh" // FormTracer oracle (refshim/)
#include "ZA_num_kernel.hh"  // numeric backend (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main() {
  using Form = DiFfRG::ZA_form_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::ZA_num_kernel;
  DressingSet d;
  auto cF = [&](double l1, double c1, double p, double k) {
    return std::real(Form::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };
  auto cN = [&](double l1, double c1, double p, double k) {
    return std::real(Num::kernel(l1, c1, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.Mq));
  };

  std::mt19937_64 rng(31337);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999);
  const int Np = 200000;
  double pw = 0;
  for (int i = 0; i < Np; ++i) {
    double l1 = U(rng), c1 = Uc(rng), p = U(rng), k = U(rng);
    double r = cF(l1, c1, p, k);
    pw = std::max(pw, std::fabs(cN(l1, c1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  auto ang1 = [&](auto kfn, double l1, double p, double k) {
    const int M = 321; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) { double c1 = -1.0 + i * h, w = (i == 0 || i == M - 1) ? 0.5 : 1.0; s += w * kfn(l1, c1, p, k); }
    return s * h;
  };
  double ie = 0;
  for (double l1 : {0.6, 1.3, 2.4}) for (double p : {0.7, 1.8}) for (double k : {0.5, 1.4}) {
    double r = ang1(cF, l1, p, k);
    ie = std::max(ie, std::fabs(ang1(cN, l1, p, k) - r) / (1e-300 + std::fabs(r)));
  }
  std::printf("ZA numeric vs FORM:  pointwise=%.3e  cos1-integrated=%.3e (physical check)\n", pw, ie);
  bool ok = ie < 1e-8;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
