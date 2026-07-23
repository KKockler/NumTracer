// Validate the NUMERIC backend for the PURE-GAUGE ZA4 four-gluon vertex flow vs the FORM oracle,
// PER TOPOLOGY. 3 loop angles.
//
// Why per topology: at full strength this flow is dominated by its ZA3^4 box, which is correct.
// The ZA3^2.ZA4 triangle and ZA4^2 bubble are the terms that were measured wrong (~0.55
// scale-relative) in NUMTRACER-ZA4-vertex-bug.md, and a single full-strength compare hides them.
// That is exactly why the pre-existing compare_za4_num.cpp (QCD quark box) passed throughout.
//
// The decomposition. Scaling the ZA3 dressing by s3 and ZA4 by s4, the integrand is exactly
//     val(s3,s4) = G + B*s3^4 + T*s3^2*s4 + U*s4^2
// (G = the ghost box, which carries neither dressing), so four evaluations separate all four:
//     G = val(0,0)                    ghost box
//     B = val(1,0) - G                ZA3^4 box
//     U = val(0,1) - G                ZA4^2 bubble      <-- was 0.54
//     T = val(1,1) - val(1,0) - val(0,1) + G   ZA3^2.ZA4 triangles  <-- was 0.55
//
// Grading: SCALE-relative (max|a-b| / max|b| over the angle grid), never pointwise-relative --
// the triangle crosses zero, where a relative error is meaningless (its FORM/NT ratio swings
// -379..+1774 while the scale-relative error is a flat 0.55). Keep the box as the pointwise
// control: it shares projector, frame, angles and routing with the broken terms, so if routing
// ever changed under us the box would break first.
#include "YM_ZA4_kernel.hh"      // FormTracer oracle (refshim/), from gen_ym_za4_form.wls
#include "YM_ZA4_num_kernel.hh"  // numeric backend   (gen/),     from gen_ym_za4_numeric.wls
#include "shim.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

int main() {
  using Form = DiFfRG::YM_ZA4_form_kernel<DiFfRG::ShimRegulator>;
  using Num = DiFfRG::YM_ZA4_num_kernel;

  // s3/s4 in {0,1} by swapping the ZA3 / ZA4 dressing for an identically-zero one.
  auto mk = [](bool s3, bool s4) {
    DressingSet d;
    if (!s3) d.ZA3 = DiFfRG::Fn{&shimdress::zero};
    if (!s4) d.ZA4 = DiFfRG::Fn{&shimdress::zero};
    return d;
  };
  auto cForm = [](const DressingSet &d, double l1, double c1, double c2, double ph, double p, double k) {
    return Form::kernel(l1, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.dtZc, d.Zc, d.dtZA, d.ZA);
  };
  auto cNum = [](const DressingSet &d, double l1, double c1, double c2, double ph, double p, double k) {
    return std::real(Num::kernel(l1, c1, c2, ph, p, k, d.ZA3, d.ZAcbc, d.ZA4, d.dtZc, d.Zc, d.dtZA, d.ZA));
  };

  // 3-angle integration on a fixed grid (trapezoid in cos1/cos2, midpoint in phi).
  auto ang3 = [&](auto fn, const DressingSet &d, double l1, double p, double k) {
    const int M = 17, Mp = 24;
    double s = 0, h = 2.0 / (M - 1), hp = 6.283185307179586 / Mp;
    for (int i = 0; i < M; ++i)
      for (int j = 0; j < M; ++j)
        for (int m = 0; m < Mp; ++m) {
          double c1 = -1 + i * h, c2 = -1 + j * h, ph = m * hp;
          double w = ((i == 0 || i == M - 1) ? 0.5 : 1.0) * ((j == 0 || j == M - 1) ? 0.5 : 1.0);
          s += w * fn(d, l1, c1, c2, ph, p, k);
        }
    return s * h * h * hp;
  };

  const DressingSet d00 = mk(false, false), d10 = mk(true, false), d01 = mk(false, true), d11 = mk(true, true);
  struct Topo {
    const char *name;
    double num, form;
  };

  // scale-relative error per topology, accumulated over a small (l1,p,k) table
  double maxAbs[4] = {}, maxScale[4] = {};
  const char *names[4] = {"ghost box (G)", "ZA3^4 box (B)", "ZA4^2 bubble (U)", "ZA3^2.ZA4 triangle (T)"};

  for (double l1 : {0.6, 1.3})
    for (double p : {0.7, 1.8})
      for (double k : {0.5, 1.4}) {
        auto ev = [&](auto fn) {
          double v00 = ang3(fn, d00, l1, p, k), v10 = ang3(fn, d10, l1, p, k);
          double v01 = ang3(fn, d01, l1, p, k), v11 = ang3(fn, d11, l1, p, k);
          return std::array<double, 4>{v00, v10 - v00, v01 - v00, v11 - v10 - v01 + v00};
        };
        auto N = ev(cNum), F = ev(cForm);
        for (int t = 0; t < 4; ++t) {
          maxAbs[t] = std::max(maxAbs[t], std::fabs(N[t] - F[t]));
          maxScale[t] = std::max(maxScale[t], std::fabs(F[t]));
        }
      }

  std::printf("pure-YM ZA4: numeric vs FORM, per topology (scale-relative, 3-angle-integrated)\n");
  bool ok = true;
  for (int t = 0; t < 4; ++t) {
    double e = maxAbs[t] / (1e-300 + maxScale[t]);
    ok = ok && (e < 1e-8);
    std::printf("  %-24s  scale-rel = %.3e   (scale %.3e)%s\n", names[t], e, maxScale[t], e < 1e-8 ? "" : "   <== FAIL");
  }

  // Pointwise control on the FULL flow: the box dominates it, so this is the routing check --
  // if routing ever diverged between the two kernels this breaks even when the topologies agree.
  double pw = 0, pwScale = 0;
  for (double l1 : {0.6, 1.3})
    for (double c1 : {-0.7, 0.3})
      for (double c2 : {-0.2, 0.8})
        for (double ph : {0.4, 3.1})
          for (double p : {0.7, 1.8})
            for (double k : {0.5, 1.4}) {
              double f = cForm(d11, l1, c1, c2, ph, p, k);
              pw = std::max(pw, std::fabs(cNum(d11, l1, c1, c2, ph, p, k) - f));
              pwScale = std::max(pwScale, std::fabs(f));
            }
  double pwe = pw / (1e-300 + pwScale);
  ok = ok && (pwe < 1e-8);
  std::printf("  %-24s  scale-rel = %.3e%s\n", "pointwise (routing)", pwe, pwe < 1e-8 ? "" : "   <== FAIL");

  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
