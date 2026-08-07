// step-05 — Hand-coding a full diagram, end to end.
//
// This is the capstone: the quark self-energy's Dirac-trace numerator,
//   T_num = tr[ p/ gamma^mu q/ gamma^nu ] P_{mu nu}(l),   q = l - p,
// assembled by hand with the numeric engine and validated two ways at one
// frame point:
//   (1) the numeric path  — describe the chain + Lorentz net as tokens,
//   contract numerically; (2) the closed form   — 4 p ( -3 c l1 + p
//   + 2 c^2 p ), the analytic value for q = l - p.
// The two agree, which is what pins the engine's algebra down.
#include <numtracer.hpp> // the whole NumTracer API — here: GlobalEnv/GenProg/emit_cpp/emit_env_layout + numeric_value/to_genprog/nprojT/dslash/dgamma

#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

namespace nm = numtracer::numeric;
namespace net = numtracer::network;
using numtracer::Cx;

// Name the axis labels. A single unscoped enum keeps every label distinct
// (auto-numbered). rho and sigma are the Lorentz dummies summed into the p- and
// q-slashes; mu and nu are the free gamma legs the projector P_{mu nu}(l) ties
// together.
enum {
  rho,
  mu,
  sigma,
  nu // Lorentz indices
};

int main() {
  // One-angle frame: p along axis 0, l at angle theta in the 0-1 plane, q = l -
  // p.
  // @snip begin: frame
  const double p = 1.3, l0 = 0.5, l1y = 0.7; // |p|, and l = (l0, l1y, 0, 0)
  const double pvec[4] = {p, 0, 0, 0};
  const double lvec[4] = {l0, l1y, 0, 0};
  const double qvec[4] = {l0 - p, l1y, 0, 0}; // q = l - p
  const double l2 = l0 * l0 + l1y * l1y;      // l^2
  // @snip end: frame

  // ---- (1) the numeric path
  // ---------------------------------------------------------------- Momenta
  // are plain numbers here, so the components are CONSTANT polynomials (nsym =
  // 0 symbols). comp[vid][m] is component m of momentum vid:  p -> 0, q -> 1, l
  // -> 2.
  const int nsym = 0;
  nm::LorentzEnv env(nsym);
  std::vector<std::array<nm::MPoly, 4>> comp(3);
  auto setv = [&](int vid, const double v[4]) {
    for (int m = 0; m < 4; ++m)
      comp[static_cast<std::size_t>(vid)][static_cast<std::size_t>(m)] =
          env.constant(Cx{v[m], 0});
  };
  setv(0, pvec);
  setv(1, qvec);
  setv(2, lvec);

  // @snip begin: tokens
  // The closed chain p/ gamma^mu q/ gamma^nu, with the two free gammas on
  // Lorentz legs mu, nu.
  net::DiracNet chain = {net::dslash({{1.0, 0}}), net::dgamma(mu),
                         net::dslash({{1.0, 1}}), net::dgamma(nu)};
  // The Lorentz network is the transverse projector P_{mu nu}(l): legs mu/nu,
  // momentum vid 2, and its 1/l^2 factor is tracked as inverse-atom id 0 (its
  // value supplied in atomDen below).
  nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nprojT(mu, nu, {{1.0, 2}}, 0)}}};
  nm::MPoly l2poly = env.zero();
  for (int m = 0; m < 4; ++m)
    l2poly = l2poly + comp[2][static_cast<std::size_t>(m)] *
                          comp[2][static_cast<std::size_t>(m)];
  std::vector<nm::MPoly> atomDen = {
      l2poly}; // atom 0 = l^2 (so its reciprocal is 1/l^2)

  nm::MPoly tr = env.numeric_value(chain, lor, comp, atomDen);
  Cx num = nm::eval(tr, /*symbols*/ {}, /*atom values*/ {1.0 / l2});
  const double numeric_val = num.re;
  // @snip end: tokens

  // ---- (2) the closed form
  // -----------------------------------------------------------------
  const double l1 = std::sqrt(l2), c = l0 / l1;
  const double closed = 4.0 * p * (-3.0 * c * l1 + p + 2.0 * c * c * p);

  std::printf("numeric  T_num = %.6f\n", numeric_val);
  std::printf(
      "closed   T_num = %.6f   (4 p(-3 c l1 + p + 2 c^2 p), q = l - p)\n",
      closed);

  // ---- (3) lower the SAME contraction to an optimized C++ kernel
  // -----------------------------------------------------------------
  // Above, the momenta were plain numbers, so the trace is a CONSTANT MPoly and would lower to a
  // bare number. To get a real *kernel* we redo the contraction with the momentum COMPONENTS as
  // symbols (ns = 3: 0 -> |p|, 1 -> l0, 2 -> l1y); the same DiracNet/Lorentz-net tokens (they refer
  // to momentum ids, not values) then contract to a polynomial in those symbols. `to_genprog`
  // Horner-factors + CSEs it into a straight-line real program over a shared symbol env `g`, and
  // `emit_cpp` prints the exact `double T_num(const double* f)` a generated kernel would carry.
  // @snip begin: lower
  const int ns = 3;
  nm::LorentzEnv envs(ns);
  auto V = [&](int i) { return envs.var(i); };            // the i-th symbol
  auto K = [&](double v) { return envs.constant(Cx{v, 0}); };
  std::vector<std::array<nm::MPoly, 4>> sc(3);
  sc[0] = {V(0), K(0), K(0), K(0)};                                 // p = (|p|, 0, 0, 0)
  sc[2] = {V(1), V(2), K(0), K(0)};                                 // l = (l0, l1y, 0, 0)
  sc[1] = {sc[2][0] - sc[0][0], sc[2][1], K(0), K(0)};              // q = l - p
  std::vector<nm::MPoly> satomDen = {sc[2][0] * sc[2][0] + sc[2][1] * sc[2][1]}; // atom 0 = l^2

  nm::MPoly trs = envs.numeric_value(chain, lor, sc, satomDen);
  net::GlobalEnv g;
  net::GenProg prog = nm::to_genprog(trs, g);
  std::printf("\nlowered kernel (symbolic components |p|, l0, l1y):\n");
  net::emit_env_layout(std::cout, g);   // which f[i] holds which symbol / 1/l^2
  net::emit_cpp(std::cout, prog, "T_num"); // -> double T_num(const double* f) { ... }
  // @snip end: lower

  // The lowered polynomial is directly runnable in-process too: eval it at the tutorial's frame
  // point and check it still equals the closed form (this keeps section (3) gated by the test).
  const double lowered_val = nm::eval(trs, {p, l0, l1y}, {1.0 / l2}).re;
  std::printf("\nsymbolic T_num = %.6f   (eval of the lowered polynomial at the same point)\n",
              lowered_val);

  const bool ok = std::fabs(numeric_val - closed) < 1e-10 &&
                  std::fabs(lowered_val - closed) < 1e-10;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
