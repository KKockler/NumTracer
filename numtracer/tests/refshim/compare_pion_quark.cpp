// Validate the NUMERIC matrix-product kernel for the Phase-2 pion<-quark-loop->pion Yukawa flow
// (Pion_quark) against the independent FormTracer (FORM) oracle — the two tensor-trace engines
// computing the SAME one-loop pion self-energy. The pion self-energy from a closed massless quark
// line with the Yukawa coupling h*qbar*i*gamma5*(tau*pi)*q and a fermion regulator R_F ~ q-slash,
// traced by NumTracer in one shot into a Dirac trace, a FUNDAMENTAL-flavour trace tr(tau^a tau^b)
// (SU(2)) and a FUNDAMENTAL-colour trace tr(1)=Nc. Both kernels bake SU(3)/SU(2) and share the
// signature kernel(l1, cos1, p, k), so this is a pointwise/angular compare. As for the other vacuum
// flows the pointwise value can differ by a loop-routing term that is odd in cos1 and integrates to
// zero — the physical check is the cos1-integrated value.
#include "Pion_quark_form_kernel.hh" // FormTracer oracle (refshim/)
#include "Pion_quark_num_kernel.hh"  // numeric backend  (gen/)
#include "shim.hpp"

#include <cmath>
#include <cstdio>
#include <random>

int main()
{
  using Form = DiFfRG::Pion_quark_form_kernel<DiFfRG::ShimRegulator>;
  using Num = Pion_quark_num_kernel;

  auto cForm = [&](double l1, double c1, double p, double k) { return Form::kernel(l1, c1, p, k); };
  auto cNum = [&](double l1, double c1, double p, double k) { return std::real(Num::kernel(l1, c1, p, k)); };

  std::mt19937_64 rng(20260703);
  std::uniform_real_distribution<double> U(0.1, 4.0), Uk(0.3, 3.0), Uc(-0.99, 0.99);
  const int Np = 200000;
  double pwF = 0;
  for (int i = 0; i < Np; ++i) {
    double l1 = U(rng), c1 = Uc(rng), p = U(rng), k = Uk(rng);
    double r = cForm(l1, c1, p, k);
    double n = cNum(l1, c1, p, k);
    pwF = std::max(pwF, std::fabs(n - r) / (1e-300 + std::fabs(r)));
  }

  // cos1-integrated check (the physical comparison: odd-in-cos1 routing terms cancel)
  auto ang1 = [&](auto kfn, double l1, double p, double k) {
    const int M = 81;
    double s = 0, h = 2.0 / (M - 1);
    for (int i = 0; i < M; ++i) {
      double c1 = -1.0 + i * h, w = (i == 0 || i == M - 1) ? 0.5 : 1.0;
      s += w * kfn(l1, c1, p, k);
    }
    return s * h;
  };
  double ieF = 0;
  for (double l1 : {0.6, 1.3, 2.4})
    for (double p : {0.4, 1.7})
      for (double k : {0.5, 1.4}) {
        double r = ang1(cForm, l1, p, k);
        double g = ang1(cNum, l1, p, k);
        ieF = std::max(ieF, std::fabs(g - r) / (1e-300 + std::fabs(r)));
      }

  std::printf("[pion <- quark loop -> pion (Yukawa) — numeric vs FORM]\n");
  std::printf("  pointwise=%.3e  cos1-integrated=%.3e (physical check)\n", pwF, ieF);
  bool ok = ieF < 1e-8;
  std::printf("%s\n", ok ? "ALL TESTS PASSED" : "TESTS FAILED");
  return ok ? 0 : 1;
}
