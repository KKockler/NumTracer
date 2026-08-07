// step-01 — hello, tensor network.
//
// The smallest possible use of the engine, with NO physics vocabulary: contract a couple of
// index networks over 4-vectors and read off the scalar they collapse to. Everything here is
// pure C++ against the header-only library — no Mathematica, no external tools.
//
// The one rule the engine runs on: two indices are summed together exactly when they carry the
// same integer label (Einstein convention). A network with no free labels left is a scalar.
//
// We do two contractions and check each against its closed form:
//   (1) a . b        = a_mu delta^{mu nu} b_nu      (a metric tying two vectors)
//   (2) a . P(k) . a = a^2 - (a.k)^2 / k^2          (a transverse projector on k)
#include <numtracer.hpp> // the whole NumTracer API — here: nvec / nmet / nprojT (NNet) + numeric_value

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
namespace net = numtracer::network;
using numtracer::Cx;

// Name the index labels with one unscoped enum: every label is auto-numbered, hence distinct, so
// two factors contract only where we deliberately reuse a label. The vector ids (a = 0, b = 1,
// k = 2) are separate plain integers that index the component table below.
enum { mu, nu }; // the two indices we will sum over

// The three 4-vectors, as plain numbers. We reuse them both to build the network and to compute
// the closed-form check, so the two can never silently drift apart.
static constexpr double a[4] = {1.0, 0.5, -0.3, 0.2};
static constexpr double b[4] = {0.8, -0.4, 1.2, 0.1};
static constexpr double k[4] = {0.5, 0.7, 0.0, 0.0};

int main() {
  // A "frame" fixes every vector's four components. Here they are literally numbers, so each
  // component is a CONSTANT polynomial (nsym = 0: an empty symbol space). comp[vid][c] is
  // component c (0..3) of vector vid.
  // @snip begin: frame
  const int nsym = 0;
  nm::LorentzEnv env(nsym);
  auto vec4 = [&env](const double v[4]) {
    return std::array<nm::MPoly, 4>{env.constant({v[0], 0}), env.constant({v[1], 0}),
                                    env.constant({v[2], 0}), env.constant({v[3], 0})};
  };
  std::vector<std::array<nm::MPoly, 4>> comp = {vec4(a), vec4(b), vec4(k)}; // vids 0, 1, 2
  // @snip end: frame

  // @snip begin: dot
  // --- (1) a . b : a metric delta_{mu nu} tying vector a (leg mu) to vector b (leg nu) -------------
  // A Lorentz network (NNet) is a sum of terms; here one term, coefficient 1, with three factors.
  // Sharing mu between a and the metric, and nu between the metric and b, sums both indices away.
  nm::NNet dot = {nm::NTerm{Cx{1, 0}, {nm::nvec(mu, {{1.0, 0}}),    // a on index mu
                                        nm::nmet(mu, nu),            // delta_{mu nu}
                                        nm::nvec(nu, {{1.0, 1}})}}};  // b on index nu
  // Empty Dirac chain (this network is pure-Lorentz). The result is one MPoly; with numeric
  // components it is a constant, and eval reads the number off (no symbols, no inverse atoms).
  const double ab = nm::eval(env.numeric_value(net::DiracNet{}, dot, comp, {}), {}, {}).re;
  // @snip end: dot

  double abClosed = 0;
  for (int c = 0; c < 4; ++c) abClosed += a[c] * b[c];

  // @snip begin: proj
  // --- (2) a . P(k) . a : the transverse projector P(k)_{mu nu} = delta_{mu nu} - k_mu k_nu/k^2 ----
  // P(k) carries its 1/k^2 as inverse-"atom" id 0; the engine needs that atom's denominator k^2.
  nm::MPoly k2 = env.constant({0, 0});
  for (int c = 0; c < 4; ++c) k2 = k2 + comp[2][c] * comp[2][c];
  nm::NNet proj = {nm::NTerm{Cx{1, 0}, {nm::nvec(mu, {{1.0, 0}}),          // a on mu
                                         nm::nprojT(mu, nu, {{1.0, 2}}, 0),  // P(k)_{mu nu}
                                         nm::nvec(nu, {{1.0, 0}})}}};        // a on nu

  double k2v = 0, a2 = 0, ak = 0;
  for (int c = 0; c < 4; ++c) k2v += k[c] * k[c], a2 += a[c] * a[c], ak += a[c] * k[c];
  const double apa = nm::eval(env.numeric_value(net::DiracNet{}, proj, comp, {k2}),
                              /*symbols*/ {}, /*atom values*/ {1.0 / k2v})
                         .re;
  // @snip end: proj
  const double apaClosed = a2 - ak * ak / k2v;

  std::printf("a . b        = %g   (expect %g)\n", ab, abClosed);
  std::printf("a . P(k) . a = %g   (expect a^2 - (a.k)^2/k^2 = %g)\n", apa, apaClosed);

  const bool ok = std::fabs(ab - abClosed) < 1e-12 && std::fabs(apa - apaClosed) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
