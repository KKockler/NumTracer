// Tutorial 6 — Per-flavour and per-component dressings: the group-diagonal fold.
//
// The flows so far dress every quark the same way (one propagator dressing for the whole
// colour/flavour multiplet). Physics often needs more: the u and d quarks dressed differently
// (broken isospin), or a single colour direction dressed on its own (a gluon condensate).
//
// The mechanism is the SU(N) colour/flavour fold. A closed SU(N) loop normally collapses to one
// flavour-BLIND number (its dimension N) — which can carry only ONE dressing. Making the loop's
// delta a *group-diagonal dressing* — diag(D_1, ..., D_N) instead of a plain delta — folds the
// same loop to a POLYNOMIAL  Sum_a D_a  over independently-named runtime dressings, WITHOUT
// splitting the trace into one diagram per component. The Dirac/Lorentz trace multiplying it is
// still computed once.
//
// Two entry points (network/sun_net.hpp):
//   sun_value_cx(net)      -- plain fold: a fully-contracted net -> one number (Cx).
//   sun_value_dressed(net) -- diagonal fold: a net with a diagFund/diagAdj factor -> a SUNPoly,
//                             i.e. Sum_t coeff_t * Prod D^dr. Each term is SUNTerm{coeff, dress},
//                             where `dress` lists the dressing-ids in that monomial. A net with no
//                             diagonal factor comes back as a single constant term == sun_value_cx.
//
// A diag factor is a delta tagged with a comp2dr map: comp2dr[v] is the dressing-id for component
// v, or -1 to DROP that component (it contributes nothing — no dead terms). Components are 0-based
// here (1..N fundamental, 1..N^2-1 adjoint in the physics, minus one for the index).
//
// We validate the two invariants the physics rests on: the COLLAPSE (all dressings equal recovers
// the flavour-blind number the plain fold gives) and the DROP (a -1 component == dressing it with
// a zero function).
#include <numtracer.hpp> // the whole NumTracer API — here: sun_value_dressed + SUN::diagFund / diagAdj / deltaFund / deltaAdj

#include <cmath>
#include <cstdio>
#include <vector>

using namespace numtracer;          // Cx, approx
namespace net = numtracer::network; // SUNPoly / SUNTerm / sun_value_dressed / SUN — as in the other tutorials

// Name the loop's SU(N) index labels. A single unscoped enum keeps every label distinct
// (auto-numbered), so axes contract iff their labels are equal.
enum { i, j }; // the two ends of the closed delta loop

// Evaluate a SUNPoly at a dressing-id -> value assignment D:  Sum_t coeff_t * Prod_{id in dress} D(id).
static Cx eval_poly(const net::SUNPoly &p, double (*D)(int)) {
  Cx s{0, 0};
  for (const net::SUNTerm &t : p) {
    Cx c = t.coeff;
    for (int id : t.dress) c = c * Cx{D(id), 0.0};
    s = s + c;
  }
  return s;
}

int main() {
  bool ok = true;

  // ---- A. Fundamental: a u/d isospin doublet (SU(2) flavour) --------------------------------
  //
  // A closed flavour delta-loop over the doublet, its delta made flavour-DIAGONAL: component 0 (u)
  // carries dressing-id 0, component 1 (d) carries id 1. diagFund(N, i, j, comp2dr) is that delta;
  // deltaFund(N, j, i) closes the loop. The fold gives the SUNPoly  D_u + D_d.
  const net::SUNPoly ud = net::sun_value_dressed(
      {net::SUN::diagFund(2, i, j, {0, 1}), net::SUN::deltaFund(2, j, i)});

  auto broken = [](int id) { return id == 0 ? 2.0 : 5.0; }; // D_u = 2, D_d = 5 (broken isospin)
  auto ones = [](int) { return 1.0; };                      // all dressings equal to 1

  const Cx ud_broken = eval_poly(ud, broken); // 2 + 5 = 7
  const Cx ud_blind = eval_poly(ud, ones);    // 1 + 1 = 2 = N_f  (the COLLAPSE)

  // DROP: the same loop dressing ONLY component 0 (id -1 drops component 1) folds to just D_u.
  const net::SUNPoly u_only = net::sun_value_dressed(
      {net::SUN::diagFund(2, i, j, {0, -1}), net::SUN::deltaFund(2, j, i)});
  const Cx u_drop = eval_poly(u_only, broken); // D_u = 2

  std::printf("A. fundamental u/d doublet (SU(2) flavour)\n");
  std::printf("   D_u + D_d           = %g   (D_u=2, D_d=5 -> 7)\n", ud_broken.re);
  std::printf("   collapse D_u=D_d=1  = %g   (-> flavour-blind N_f = 2)\n", ud_blind.re);
  std::printf("   drop d (only u)     = %g   (-> D_u = 2)\n", u_drop.re);
  ok = ok && approx(ud_broken, Cx{7, 0}) && approx(ud_blind, Cx{2, 0}) && approx(u_drop, Cx{2, 0});

  // ---- B. Adjoint: a gluon condensate on the Cartan directions (SU(3) colour) ---------------
  //
  // The adjoint has N^2-1 = 8 components. diagAdj(N, i, j, comp2dr) dresses them individually. A
  // condensate lives in the CARTAN directions -- lambda_3, lambda_8 for SU(3), i.e. components 3
  // and 8 (0-based indices 2 and 7). We dress those two and DROP the other six.
  std::vector<int> cartan(8, -1); // start with everything dropped
  cartan[2] = 0;                  // lambda_3 -> dressing-id 0  (Z_3)
  cartan[7] = 1;                  // lambda_8 -> dressing-id 1  (Z_8)
  const net::SUNPoly cond = net::sun_value_dressed(
      {net::SUN::diagAdj(3, i, j, cartan), net::SUN::deltaAdj(3, j, i)});
  const Cx cond_blind = eval_poly(cond, ones); // Z_3 + Z_8 at 1 -> 2

  // The colour-BLIND adjoint loop, every component dressed by one id, is the collapse baseline:
  // all-equal -> N^2-1 = 8 (what the plain sun_value_cx delta-loop would give).
  std::vector<int> all8(8);
  for (int a = 0; a < 8; ++a) all8[a] = a;
  const net::SUNPoly full = net::sun_value_dressed(
      {net::SUN::diagAdj(3, i, j, all8), net::SUN::deltaAdj(3, j, i)});
  const Cx full_blind = eval_poly(full, ones); // 8

  // DROP == zero-dressing: the Cartan-only poly equals the full poly with the other six dressings
  // set to zero. Here: full evaluated with D(id)=1 only for the Cartan ids {2,7}, else 0 -> 2.
  auto cartan_mask = [](int id) { return (id == 2 || id == 7) ? 1.0 : 0.0; };
  const Cx full_masked = eval_poly(full, cartan_mask); // 2, matching cond_blind

  std::printf("B. adjoint gluon condensate on the Cartan (SU(3) colour)\n");
  std::printf("   full loop, all Z=1  = %g   (-> N^2-1 = 8)\n", full_blind.re);
  std::printf("   Cartan {3,8}, Z=1   = %g   (only lambda_3, lambda_8 survive -> 2)\n", cond_blind.re);
  std::printf("   full with rest=0    = %g   (drop == zero-dressing -> matches Cartan)\n", full_masked.re);
  ok = ok && approx(full_blind, Cx{8, 0}) && approx(cond_blind, Cx{2, 0}) &&
       approx(full_masked, cond_blind);

  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
