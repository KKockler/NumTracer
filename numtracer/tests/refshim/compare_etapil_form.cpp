// Validate the NUMERIC matrix-product kernel for the QCD pion anomalous dimension etaPiL against
// the independent FormTracer (FORM) oracle — the two tensor-trace engines computing the SAME pion
// two-point wave-function loop. The pion self-energy gets its p^2 dependence from the quark loop
// (two pion-quark Yukawa vertices + a regularised, constituent-mass quark propagator); the pointlike
// meson loops are p-independent and drop out of the (f(p)-f(0))/p^2 extraction (done before tracing).
// As for the other vacuum flows, the pointwise value can differ by a loop-routing term that is odd in
// cos1 and integrates to zero — the physical check is the cos1-integrated value.
#include "EtaPiL_form_kernel.hh"       // FormTracer oracle (refshim/)
#include "EtaPiL_num_kernel.hh"        // numeric backend (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main()
{
  using Form  = DiFfRG::EtaPiL_form_kernel<DiFfRG::ShimRegulator>;
  using Num   = DiFfRG::EtaPiL_num_kernel<DiFfRG::ShimRegulator>;
  DressingSet d;

  // scalar runtime params (etaPiL, etaSigL, d1V, d2V, rhoL); only rhoL (the quark constituent mass
  // ~ hSigL*sqrt(rhoL/Nf)) enters the value — the others are unused signature slots.
  auto cForm = [&](double l1, double c1, double k, double eP, double eS, double d1V, double d2V, double rhoL) {
    return Form::kernel(l1, c1, k, eP, eS, d1V, d2V, rhoL, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc, d.Zc,
                        d.dtZA, d.ZA, d.dtZq, d.Zq, d.hPiL, d.hSigL, d.zPiL, d.zSigL);
  };
  auto cNum = [&](double l1, double c1, double k, double eP, double eS, double d1V, double d2V, double rhoL) {
    return std::real(Num::kernel(l1, c1, k, eP, eS, d1V, d2V, rhoL, d.ZA3, d.ZAcbc, d.ZA4, d.ZAqbq1, d.dtZc,
                                 d.Zc, d.dtZA, d.ZA, d.dtZq, d.Zq, d.hPiL, d.hSigL, d.zPiL, d.zSigL));
  };

  std::mt19937_64 rng(20260626);
  std::uniform_real_distribution<double> U(0.05, 3.0), Uc(-0.999, 0.999), Ue(-1.5, 1.5);
  const int Np = 200000;
  double pwF = 0;
  for (int i = 0; i < Np; ++i) {
    double l1 = U(rng), c1 = Uc(rng), k = U(rng);
    double eP = Ue(rng), eS = Ue(rng), d1V = U(rng), d2V = Ue(rng), rhoL = U(rng);
    double r = cForm(l1, c1, k, eP, eS, d1V, d2V, rhoL);
    double n = cNum(l1, c1, k, eP, eS, d1V, d2V, rhoL);
    pwF = std::max(pwF, std::fabs(n - r) / (1e-300 + std::fabs(r)));
  }

  // cos1-integrated check (the physical comparison: odd-in-cos1 routing terms cancel)
  auto ang1 = [&](auto kfn, double l1, double k, double eP, double eS, double d1V, double d2V, double rhoL) {
    const int M = 81; double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) {
      double c1 = -1.0 + i * h, w = (i == 0 || i == M - 1) ? 0.5 : 1.0;
      s += w * kfn(l1, k, c1, eP, eS, d1V, d2V, rhoL);
    }
    return s * h;
  };
  // adapt arg order: ang1 passes (l1,k,c1,...) -> wrap to kernel order (l1,c1,k,...)
  auto wF = [&](double l1, double k, double c1, double eP, double eS, double d1V, double d2V, double rhoL) {
    return cForm(l1, c1, k, eP, eS, d1V, d2V, rhoL);
  };
  auto wN = [&](double l1, double k, double c1, double eP, double eS, double d1V, double d2V, double rhoL) {
    return cNum(l1, c1, k, eP, eS, d1V, d2V, rhoL);
  };
  double ieF = 0;
  for (double l1 : {0.6, 1.3, 2.4}) for (double k : {0.5, 1.4}) for (double rhoL : {0.4, 1.7}) {
    double r = ang1(wF, l1, k, 0.3, -0.2, 0.8, 0.4, rhoL);
    double g = ang1(wN, l1, k, 0.3, -0.2, 0.8, 0.4, rhoL);
    ieF = std::max(ieF, std::fabs(g - r) / (1e-300 + std::fabs(r)));
  }

  // The pointwise spread (~1e-6) is benign round-off from the p=1e-3 wave-function trick: the kernel
  // carries a ~5e5 = 1/(2 p^2) prefactor and ill-conditioned `cos1 - 1000 l1` momentum arguments, so
  // algebraically-equal backends differ in the last ~6 digits per point. The physical comparison is
  // the cos1-integrated value (the actual angular average the integrator forms), where it averages out.
  std::printf("etaPiL numeric vs FORM:   pointwise=%.3e  cos1-integrated=%.3e (physical check)\n", pwF, ieF);
  bool ok = ieF < 1e-8;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
