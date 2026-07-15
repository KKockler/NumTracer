// Tutorial 7 — Dressed propagator numerators: one trace instead of 2^D diagrams.
//
// A dressed quark propagator numerator is a SUM of Dirac structures whose coefficients are
// *runtime* dressings, e.g.
//
//     S(p) = Mq · 𝟙  +  Z(p) · p̸        (a mass term + a wave-function-dressed slash)
//
// A diagram with D such dressed numerators is, if you distribute the sum, 2^D separate traces —
// one per choice of structure in each numerator. NumTracer does NOT distribute. It keeps each
// numerator EAGER (a "slot" = the list of its structure options) and contracts the chain ONCE,
// collecting the result into a DPoly: a polynomial whose *variables* are the dressing calls and
// whose *coefficients* are the kinematic MPolys the engine already computes (mpoly.hpp). The
// dressings ride along as opaque atom-ids and never enter the trace arithmetic — so the Dirac /
// Lorentz work is done a single time no matter how many structures each numerator carries.
//
// This is the Dirac-side analogue of tutorial 6's colour/flavour fold (there a group-diagonal δ
// folds to a SUNPoly; here a dressed numerator folds to a DPoly). See dpoly.hpp for the type and
// internals/numeric-engine.md for how DPoly wraps MPoly.
//
// The example: a quark bubble  tr( γ^μ · S(p) · γ^ν · S(q) )  closed by the gluon metric δ_{μν},
// with two dressed numerators S(p), S(q). Distributed, that is 2×2 = 4 traces; collected, it is
// ONE DPoly of (at most) 4 dressing monomials. We validate that the collected DPoly, evaluated at
// random kinematics and random dressing values, equals the explicit distributed sum to 1e-10.
#include <numtracer.hpp> // the whole API — here: numeric_value_dressed / DSlot / DPoly / eval + network::dgamma / dslash

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using numtracer::Cx;                 // the complex scalar
namespace nm = numtracer::numeric;   // MPoly / DPoly / DSlot / numeric_value / numeric_value_dressed / eval
namespace net = numtracer::network;  // DiracNet / dgamma / dslash / NetVal builders

static double cdiff(Cx a, Cx b) { return std::abs(a.re - b.re) + std::abs(a.im - b.im); }

int main() {
  // Frame symbols: two momenta p (symbols 0..3) and q (symbols 4..7). comp[vid] is the 4-vector of
  // MPoly components of fundamental momentum `vid` — here each component is just a bare symbol.
  const int nsym = 8; // p:0..3, q:4..7
  std::vector<std::array<nm::MPoly, 4>> comp(2);
  for (int mu = 0; mu < 4; ++mu) {
    comp[0][mu] = nm::MPoly::var(nsym, mu);     // p_mu
    comp[1][mu] = nm::MPoly::var(nsym, 4 + mu); // q_mu
  }

  // The Lorentz half: the two free gluon legs (ids 100, 101) meet through one metric δ_{100,101}.
  nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nmet(100, 101)}}};

  // The dressed numerators. A DSlot is the list of a numerator's structure options; a DSlotOpt is
  //   { coeff , dressing-atom ids , slash? , slash-momentum } .
  //   S(p): option 0 = identity 𝟙 dressed by atom 0 ("Mq");   option 1 = p̸ dressed by atom 1 ("Z(p)").
  //   S(q): option 0 = identity 𝟙 dressed by atom 0 ("Mq");   option 1 = q̸ dressed by atom 2 ("Z(q)").
  // (Atom 0 is shared: both mass terms are the SAME runtime Mq. Atoms 1,2 are the per-leg Z's.)
  nm::DSlot sP = {nm::DSlotOpt{Cx{1, 0}, {0}, /*slash*/ false, {}},
                  nm::DSlotOpt{Cx{1, 0}, {1}, /*slash*/ true, {{1.0, 0}}}}; // p̸ = 1·comp[0]
  nm::DSlot sQ = {nm::DSlotOpt{Cx{1, 0}, {0}, /*slash*/ false, {}},
                  nm::DSlotOpt{Cx{1, 0}, {2}, /*slash*/ true, {{1.0, 1}}}}; // q̸ = 1·comp[1]

  // The dressed Dirac chain, in trace order: a token is either a FIXED factor (dtfix) or a SLOT
  // reference (dtslot i -> the i-th entry of the slot list below).  γ^100 · S(p) · γ^101 · S(q).
  std::vector<nm::DChainTok> dchain = {nm::dtfix(net::dgamma(100)), nm::dtslot(0),
                                       nm::dtfix(net::dgamma(101)), nm::dtslot(1)};

  // Collect: ONE contraction, no 2^D blowup. dp is the DPoly.
  nm::DPoly dp = nm::numeric_value_dressed(nsym, dchain, {sP, sQ}, lor, comp, /*atomDen*/ {});

  // Distributed reference: enumerate the 2×2 structure choices, contract each concrete (undressed)
  // chain with the ORDINARY numeric_value, and weight by the product of that choice's dressings.
  auto refTrace = [&](int cp, int cq, const std::vector<double> &x) {
    net::DiracNet c = {net::dgamma(100)};
    if (cp == 1) c.push_back(net::dslash({{1.0, 0}})); // p̸ (else identity: nothing to push)
    c.push_back(net::dgamma(101));
    if (cq == 1) c.push_back(net::dslash({{1.0, 1}})); // q̸
    return nm::eval(nm::numeric_value(nsym, c, lor, comp, {}), x, {});
  };

  std::printf("Dressed quark bubble  tr( γ^μ S(p) γ^ν S(q) ) δ_{μν},  S = Mq·𝟙 + Z·slash\n");
  std::printf("  DPoly dressing monomials (one trace collected) : %d   (distributed would be 4 traces)\n",
              dp.size());

  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  double maxerr = 0.0;
  int bad = 0;
  for (int it = 0; it < 5000; ++it) {
    std::vector<double> x(nsym);
    for (double &v : x) v = U(rng);
    // 3 entries because the slot options above reference dressing-ids {0,1,2}; nothing tells DPoly
    // "there are 3 dressings" — eval() indexes drVal BY id (drVal[id]), so it must be sized to
    // max-id + 1. drVal[0]=Mq (shared by both mass terms), drVal[1]=Z(p), drVal[2]=Z(q).
    std::vector<double> drVal = {U(rng), U(rng), U(rng)};

    // collected: evaluate the single DPoly (kinematics x, no 1/k² atoms here, dressings drVal)
    Cx collected = nm::eval(dp, x, /*atomVal*/ {}, drVal);

    // distributed: Σ over the 4 structure choices of (dressing product) · trace(choice)
    Cx dist{0, 0};
    for (int cp = 0; cp < 2; ++cp)
      for (int cq = 0; cq < 2; ++cq) {
        double w = (cp == 0 ? drVal[0] : drVal[1]) * (cq == 0 ? drVal[0] : drVal[2]);
        Cx tr = refTrace(cp, cq, x);
        dist = dist + Cx{tr.re * w, tr.im * w};
      }

    double e = cdiff(collected, dist);
    maxerr = std::max(maxerr, e);
    if (e >= 1e-10) ++bad;
  }

  const bool ok = (bad == 0);
  std::printf("  collected == distributed over 5000 random points : worst |Δ| = %.2e  (%d bad)\n", maxerr, bad);
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
