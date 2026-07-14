// The generator's sub-term dedup (numeric/trace_fold.hpp): phase A contracts each DISTINCT Dirac
// trace once, phase B folds each net's traces with its own scalars as a balanced binary-counter tree.
//
// Why this test exists. The dedup only does anything on a flow with MANY sub-terms per net, and every
// committed flow in tests/gen has exactly ONE — so the tree fold degenerates to a single push there
// and the whole suite, codegen test included, is blind to it. The fold is also the piece most likely
// to be silently wrong: a carry-propagation stack that drops or double-counts a partial still returns
// a plausible polynomial. So exercise it directly, against the naive left fold it replaces, at the
// sizes and cache regimes the real flows hit:
//
//   - fold_net == the left fold  (reassociation must not change the VALUE; MPoly/DPoly addition is
//     exact term collection over Cx, so this is an exact equality, not a tolerance)
//   - every cache regime: nCache = 0 (nothing resident, all recomputed), 1, n/2, n (all resident)
//   - a trace beyond nCache is recomputed EXACTLY ONCE (the singleton-eviction assumption: Codegen.m
//     only leaves refCount==1 traces uncached, so a recompute must never be duplicated work)
//   - both backends: MPoly (plain) and DPoly (dressed)
//   - the degenerate shapes that the binary counter gets wrong if the carry is off by one: 0, 1, 2, 3
//     and 7 terms (7 = 0b111 leaves three unmerged partials on the stack -> the leftover sweep)
#include "numtracer/numeric/trace_fold.hpp"

#include <cstdio>
#include <numeric>
#include <vector>

using namespace numtracer;
using namespace numtracer::numeric;

namespace
{
  int failures = 0;

  void check(bool ok, const char *what)
  {
    if (!ok) {
      std::printf("  FAIL: %s\n", what);
      ++failures;
    }
  }

  constexpr int NSYM = 3;

  // A distinct, non-trivial polynomial per trace id: c(k)·x0^(k%3) · x1^(k%2) + k·x2. Distinct enough
  // that a dropped or double-counted term cannot cancel out of the total by luck.
  MPoly trace_of(int k)
  {
    MPoly p = MPoly::mono(NSYM, {k % 3, k % 2, 0}, Cx{1.0 + k, 0.5 * k});
    return p + MPoly::mono(NSYM, {0, 0, 1}, Cx{static_cast<double>(k), 0.0});
  }

  bool same(const MPoly &a, const MPoly &b)
  {
    if (a.t.size() != b.t.size()) return false;
    for (std::size_t i = 0; i < a.t.size(); ++i) {
      if (!(a.t[i].first == b.t[i].first)) return false;
      if (a.t[i].second.re != b.t[i].second.re || a.t[i].second.im != b.t[i].second.im) return false;
    }
    return true;
  }

  bool same(const DPoly &a, const DPoly &b)
  {
    if (a.t.size() != b.t.size()) return false;
    for (std::size_t i = 0; i < a.t.size(); ++i) {
      if (a.t[i].first != b.t[i].first) return false;
      if (!same(a.t[i].second, b.t[i].second)) return false;
    }
    return true;
  }

  // A stand-in polynomial that records how many ELEMENTS each addition touches. `size` models the
  // term count; the cost of adding two polynomials is linear in their combined size (MPoly::operator+
  // is a merge). Summing that over the fold is the quantity the tree fold is built to keep at
  // O(n log n) — a left fold against a growing accumulator makes it O(n²).
  struct Counted {
    long size = 0;
    static long touches;
    Counted() = default;
    explicit Counted(int) {}                                  // the P(nsym) empty accumulator
    static Counted terms(long s) { Counted c; c.size = s; return c; }
  };
  long Counted::touches = 0;

  Counted operator+(const Counted &a, const Counted &b)
  {
    Counted::touches += a.size + b.size;                      // a merge touches both operands
    return Counted::terms(a.size + b.size);
  }
  Counted scale_trace(int, const Counted &p, Cx) { return p; } // found by ADL from inside fold_net

  // What fold_net replaces: Σ_j sc[j]·trace(idx[j]), left-folded against one growing accumulator.
  template <class P, class TraceFn>
  P left_fold(const std::vector<int> &idx, const std::vector<Cx> &sc, TraceFn &&trace)
  {
    P acc(NSYM);
    for (std::size_t j = 0; j < idx.size(); ++j)
      acc = acc + scale_trace(NSYM, trace(idx[j]), sc[j]);
    return acc;
  }

  // One (backend, net shape, cache regime) case.
  template <class P, class MakeTrace>
  void run_case(const char *tag, std::size_t n, long nCache, MakeTrace &&mk)
  {
    std::vector<int> idx(n);
    std::vector<Cx> sc(n);
    for (std::size_t j = 0; j < n; ++j) {
      idx[j] = static_cast<int>(j);            // every trace distinct => refCount 1 => recompute is 1x
      sc[j] = Cx{1.0 + 0.25 * j, -0.125 * j};  // non-uniform scalars: a swapped scalar shows up
    }

    // Phase A over the resident prefix, exactly as the generated main does it.
    std::vector<P> T = contract_traces<P>(NSYM, nCache, /*W=*/4, mk);
    check(T.size() == static_cast<std::size_t>(nCache < 0 ? 0 : nCache), "phase A table size");

    // Count recomputes so we can prove nothing beyond the cache is computed twice.
    std::vector<int> hits(n, 0);
    auto counting = [&](int k) {
      if (static_cast<std::size_t>(k) < hits.size()) ++hits[static_cast<std::size_t>(k)];
      return mk(k);
    };

    const P got = fold_net<P>(NSYM, idx, sc, T, nCache, counting);
    const P want = left_fold<P>(idx, sc, mk);

    char msg[160];
    std::snprintf(msg, sizeof msg, "%s n=%zu nCache=%ld: tree fold == left fold", tag, n, nCache);
    check(same(got, want), msg);

    for (std::size_t k = 0; k < n; ++k) {
      const int expect = (static_cast<long>(k) < nCache) ? 0 : 1; // resident => 0 recomputes, else 1
      if (hits[k] != expect) {
        std::snprintf(msg, sizeof msg, "%s n=%zu nCache=%ld: trace %zu recomputed %dx (want %dx)", tag, n,
                      nCache, k, hits[k], expect);
        check(false, msg);
        break;
      }
    }
  }

  template <class P, class MakeTrace> void run_backend(const char *tag, MakeTrace &&mk)
  {
    // 0/1/2/3 are the binary-counter edge cases; 7 = 0b111 leaves three partials for the leftover
    // sweep; 1000 is a heavy net, where the left fold's O(n^2) is what the tree is there to avoid.
    for (std::size_t n : {std::size_t{0}, std::size_t{1}, std::size_t{2}, std::size_t{3}, std::size_t{7},
                          std::size_t{64}, std::size_t{1000}}) {
      for (long nCache : {0L, 1L, static_cast<long>(n) / 2, static_cast<long>(n)}) {
        if (nCache > static_cast<long>(n)) continue;
        run_case<P>(tag, n, nCache, mk);
      }
    }
  }
} // namespace

int main()
{
  std::printf("trace_fold: tree fold vs left fold, all cache regimes\n");

  run_backend<MPoly>("MPoly", trace_of);
  run_backend<DPoly>("DPoly", [](int k) {
    // A dressed trace: the same kinematic polynomial under a dressing monomial that varies with k, so
    // the DPoly merge (which is keyed on the dressing monomial) is exercised rather than bypassed.
    DPoly d(NSYM);
    d.add(dmono_sorted(DMono{k % 2, k % 3}), trace_of(k));
    return d;
  });

  // A net that references the SAME trace repeatedly — the shape the dedup actually produces (a hot
  // trace is referenced 2-32x across the flow). Left fold and tree fold must still agree.
  {
    const std::vector<int> idx = {2, 2, 5, 2, 5, 0, 5, 5};
    const std::vector<Cx> sc = {{1, 0}, {-1, 0.5}, {2, 0}, {0.5, -0.25}, {1, 1}, {3, 0}, {-2, 0}, {0.25, 0}};
    std::vector<MPoly> T = contract_traces<MPoly>(NSYM, 6, 4, trace_of);
    check(same(fold_net<MPoly>(NSYM, idx, sc, T, 6, trace_of), left_fold<MPoly>(idx, sc, trace_of)),
          "repeated-trace net: tree fold == left fold");
  }

  // The tree must actually be a TREE. Every correct fold adds each term exactly once, so a botched
  // carry that merges only one level still returns the RIGHT VALUE — it just degenerates back into
  // the left fold this exists to replace. A value check cannot see that (verified: turning the carry
  // loop into a single `if` passes every check above). So bound the WORK instead: count the operand
  // sizes touched across all additions, which is what a left fold makes quadratic.
  {
    constexpr std::size_t n = 1024;
    std::vector<int> idx(n);
    std::vector<Cx> sc(n, Cx{1, 0});
    for (std::size_t j = 0; j < n; ++j)
      idx[j] = static_cast<int>(j);

    const std::vector<Counted> none;
    Counted::touches = 0;
    (void)fold_net<Counted>(NSYM, idx, sc, none, 0, [](int) { return Counted::terms(1); });
    const long tree = Counted::touches;

    Counted::touches = 0;
    (void)left_fold<Counted>(idx, sc, [](int) { return Counted::terms(1); });
    const long left = Counted::touches;

    // n·log2(n) ≈ 10k vs n²/2 ≈ 524k. The 4x margin keeps this from being brittle while still being
    // orders of magnitude away from a left fold.
    char msg[160];
    std::snprintf(msg, sizeof msg, "tree fold is O(n log n): %ld touches vs left fold %ld (n=%zu)", tree,
                  left, n);
    check(tree * 4 < left, msg);
  }

  // poly_bytes must actually track size — it is what the generator reports as the trace table's RAM
  // cost, and a constant would silently hide the one way this design can regress (peak RSS).
  check(poly_bytes(MPoly(NSYM)) == 0, "poly_bytes: empty is 0");
  check(poly_bytes(trace_of(4)) > poly_bytes(MPoly::constant(NSYM, Cx{1, 0})), "poly_bytes: grows with terms");

  if (failures) {
    std::printf("trace_fold: %d FAILURES\n", failures);
    return 1;
  }
  std::printf("trace_fold: all checks passed\n");
  return 0;
}
