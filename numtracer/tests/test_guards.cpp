// NumTracer — guards against silent-corruption traps on code paths no current flow exercises.
//
// These cover latent bugs the end-to-end FORM/dense oracles cannot catch because no committed flow
// triggers them:
//   (A1) network::add collapses a same-index vector sum into one compound-vector leaf whose weights
//        are REAL (Elem::vlc is std::pair<double,int>). A complex coefficient on such a term would
//        have its imaginary part silently dropped; the collapse must refuse it loudly instead.
//   (B3) an EMPTY colour network is the identity factor (1), not a UnionFind(-1)/degenerate path.
//   (C1) dirac_value's gamma-parity count must treat a Comm (sigma) as TWO gammas and a LoopSep as
//        NONE — and, since its Wick recursion can trace neither, must refuse both loudly. The old
//        `!= Gamma5` count called each ONE gamma, which inverted the verdict in both directions:
//        {sigma,gamma,gamma} (4 gammas, nonzero trace) was reported as an odd chain and returned a
//        structural ZERO, while {sigma,gamma} (3 gammas, vanishing trace) passed the even gate into
//        trace_rec, whose pair_factor read the Comm's empty vlc and silently collapsed it. This
//        function is the cross-validation ORACLE in test_numeric_contract.cpp, so a wrong verdict
//        here would "confirm" a wrong engine result.
#include "numtracer/network/dirac.hpp"
#include "numtracer/network/network.hpp"
#include "numtracer/network/sun_net.hpp"
#include "numtracer/sun/sun_data.hpp"

#include <cstdio>
#include <stdexcept>

using namespace numtracer;

int main() {
  int fail = 0;
  auto ok = [&](const char *name, bool pass) {
    std::printf("  %s: %s\n", pass ? "ok  " : "FAIL", name);
    if (!pass) ++fail;
  };

  // ---- A1: real-coefficient vector sums collapse; complex ones throw ------------------------
  {
    using namespace numtracer::network;
    // two momentum legs on the same Lorentz index, real weights → one compound-vector leaf.
    NetVal s = add(scale(Cx{2, 0}, vec(0, 5)), scale(Cx{-1, 0}, vec(0, 7)));
    ok("real vecsum collapses to one term", s.size() == 1 && s[0].e.size() == 1);
    ok("compound vector carries both momenta", !s.empty() && s[0].e[0].vlc.size() == 2);

    bool threw = false;
    try {
      (void)add(scale(Cx{0, 1}, vec(0, 5)), scale(Cx{0, 1}, vec(0, 7)));
    } catch (const std::exception &) {
      threw = true;
    }
    ok("complex vecsum coefficient throws (no silent im-drop)", threw);
  }

  // ---- B3: the empty colour network is the identity factor ----------------------------------
  {
    using namespace numtracer::network;
    const Cx e = sun_value_cx(SUNNet{});
    ok("sun_value_cx(empty) == 1", e.re == 1.0 && e.im == 0.0);

    const SUNPoly p = sun_value_dressed(SUNNet{});
    const bool unit = (p.size() == 1 && p[0].dress.empty() && p[0].coeff.re == 1.0 && p[0].coeff.im == 0.0);
    ok("sun_value_dressed(empty) == single unit term", unit);
  }

  // ---- C1: dirac_value refuses tokens its Wick recursion cannot trace ------------------------
  {
    using namespace numtracer::network;
    auto refuses = [&](const DiracNet &ch) {
      try {
        (void)dirac_value(ch, 900000);
        return false;
      } catch (const std::exception &) {
        return true;
      }
    };
    // physically 4 gammas (nonzero trace) — the old count said "odd" and returned a structural zero
    ok("dirac_value refuses {sigma,gamma,gamma}", refuses({dcomm(0, 1), dgamma(2), dgamma(3)}));
    // physically 3 gammas (vanishing trace) — the old count said "even" and fell into trace_rec
    ok("dirac_value refuses {sigma,gamma}", refuses({dcomm(0, 1), dgamma(2)}));
    ok("dirac_value refuses a LoopSep chain", refuses({dgamma(0), dloopsep(), dgamma(1)}));

    // plain gamma/slash chains must be completely unaffected by the guard + recount
    ok("even gamma chain still traces", !dirac_value({dgamma(0), dgamma(1)}, 900000).empty());
    ok("odd gamma chain still vanishes", dirac_value({dgamma(0), dgamma(1), dgamma(2)}, 900000).empty());
    ok("slash pair still traces",
       !dirac_value({dslash({{1.0, 0}}), dslash({{1.0, 4}})}, 900000).empty());
  }

  std::printf("\n%s (%d failure%s)\n", fail ? "TESTS FAILED" : "ALL TESTS PASSED", fail, fail == 1 ? "" : "s");
  return fail ? 1 : 0;
}
