/// @file trace_fold.hpp
/// @brief The generator's two contraction phases: contract each DISTINCT Dirac trace once (phase A),
///        then fold each net's traces with its own scalars (phase B).
///
/// The generator's cost is one trace contraction per `(net, sub-term)`, but the same
/// `(dnet, lnet, dch, dsl)` tuple recurs across nets and colour branches, so most of those
/// contractions recompute a trace that was already computed. Measured on dense flows: 30,807
/// contractions for 6,041 distinct traces (5.1x), and 246,456 for 32,784 (7.5x). Codegen.m therefore
/// emits a table of the distinct traces and, per net, the indices into it; this header contracts and
/// folds that table. Two independent wins:
///
///  - **the contraction shrinks by the redundancy factor** — each distinct trace is contracted once;
///  - **the parallel phase becomes a FLAT list of uniform work items.** Scheduling per *net* could
///    not use the machine: sub-terms per net are wildly skewed (max 2880 vs a median of 27), so the
///    single biggest net alone exceeded the ideal per-thread load and pinned utilisation at ~33%
///    however many cores were available.
///
/// Phase B reduces each net as a **balanced binary-counter tree** rather than a left fold against a
/// growing accumulator: heavy nets fold >1000 terms, and a left fold pays O(|acc|) on every one of
/// them. Reassociating perturbs the emitted kernel only in the last ulp (measured: worst relative
/// deviation ~7e-15), which is within the accepted tolerance for the generated coefficients.
///
/// Only the generator's main TU includes this (it is host-only: it spawns threads). The `-O0` net
/// builder units are compiled `-fno-exceptions -fno-rtti` and must not see it.
#pragma once

#include "numtracer/numeric/dpoly.hpp"
#include "numtracer/numeric/mpoly.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

namespace numtracer::numeric
{

  /// @brief Scale a contracted trace by its sub-term scalar. Overloaded so phase B is one template
  ///        over both backends: the plain path multiplies by a constant `MPoly`, the dressed path
  ///        scales every kinematic coefficient (@ref scaleCx).
  inline MPoly scale_trace(int nsym, const MPoly &p, Cx c) { return MPoly::constant(nsym, c) * p; }
  inline DPoly scale_trace(int, const DPoly &p, Cx c) { return scaleCx(p, c); }

  /// @brief Approximate heap footprint of a polynomial, for the trace-table RAM report. The table is
  ///        the one place this design trades memory for time, so the generator prints what it costs.
  inline std::size_t poly_bytes(const MPoly &p) { return p.t.size() * sizeof(std::pair<Mono, Cx>); }
  inline std::size_t poly_bytes(const DPoly &p)
  {
    std::size_t b = 0;
    for (const auto &[d, m] : p.t)
      b += d.size() * sizeof(int) + poly_bytes(m);
    return b;
  }

  /// @brief Run `f(i)` for `i` in `[0, n)` across `W` threads, dispatched off one flat atomic counter.
  ///        Work-stealing rather than a static split: even a uniform work list has a long tail, and
  ///        this is also what lets phase A ignore the per-net skew entirely. Returns the worker count
  ///        actually used. Falls back to fewer threads (down to the caller's own) if the system
  ///        refuses to spawn them.
  template <class F> unsigned parallel_flat(long n, unsigned W, F &&f)
  {
    if (n <= 0) return 0u;
    const unsigned nw = static_cast<unsigned>(std::min<long>(std::max(1u, W), n));
    std::atomic<long> next{0};
    auto work = [&] {
      long i;
      while ((i = next.fetch_add(1)) < n) f(i);
    };
    std::vector<std::thread> pool;
    pool.reserve(nw - 1);
    for (unsigned w = 1; w < nw; ++w) {
      try {
        pool.emplace_back(work);
      } catch (const std::system_error &) {
        break; // out of threads: the remaining work still runs, just on fewer of them
      }
    }
    work();
    for (auto &t : pool)
      t.join();
    return nw;
  }

  /// @brief PHASE A — contract traces `[0, nCache)` once each, in parallel over a flat work list.
  ///
  /// `nCache` is how many of the distinct traces are cached. Codegen.m orders the traces by
  /// *descending reference count*, so a cap keeps the traces that repay caching most; traces at index
  /// `>= nCache` are recomputed on demand in phase B (see @ref fold_net). Two reasons to cap:
  ///  - a trace referenced exactly once costs the same either way, so caching it is pure RAM for no
  ///    saving — Codegen.m's default `nCache` therefore excludes the singletons;
  ///  - `NT_GEN_MEMO_MAX` trims it further when memory is tight (the RAM lever: the dense flows are
  ///    memory-bound before they are compute-bound).
  ///
  /// @param trace `trace(k) -> P`, the contraction of distinct trace `k`. Must be pure and safe to
  ///        call concurrently (the numeric_value_* entry points take everything by const reference).
  template <class P, class TraceFn>
  std::vector<P> contract_traces(int nsym, long nCache, unsigned W, TraceFn &&trace)
  {
    std::vector<P> T(static_cast<std::size_t>(std::max(0L, nCache)), P(nsym));
    parallel_flat(nCache, W, [&](long k) { T[static_cast<std::size_t>(k)] = trace(static_cast<int>(k)); });
    return T;
  }

  /// @brief PHASE B — fold one net: `Σ_j sc[j] · trace(idx[j])`.
  ///
  /// Reduced as a **binary-counter tree**: partial sums are kept on a stack tagged with a rank (a
  /// power of two = how many leaves it covers), and a new term carries into the stack exactly as a
  /// binary increment does, merging only equal ranks. So every addition combines two operands of
  /// comparable size, and at most O(log n) partials are ever live — where a left fold against one
  /// growing accumulator pays O(|acc|) per term, which on the heavy nets (>1000 terms) is the
  /// difference between O(n log n) and O(n²) element touches. Exactly `n - 1` additions either way.
  ///
  /// Traces at `idx[j] >= nCache` are not resident and are recomputed here; by construction those are
  /// referenced once, so nothing is computed twice.
  template <class P, class TraceFn>
  P fold_net(int nsym, const std::vector<int> &idx, const std::vector<Cx> &sc, const std::vector<P> &T,
             long nCache, TraceFn &&trace)
  {
    std::vector<P> st;
    std::vector<std::size_t> rank;
    P recomputed;

    for (std::size_t j = 0; j < idx.size(); ++j) {
      const int k = idx[j];
      // Bind, never copy: a ternary over `const P&` and a prvalue would materialise a copy of the
      // cached polynomial on every use — and the whole point is that these are used 5-8x each.
      const P *src;
      if (k < static_cast<int>(nCache)) {
        src = &T[static_cast<std::size_t>(k)];
      } else {
        recomputed = trace(k);
        src = &recomputed;
      }

      P cur = scale_trace(nsym, *src, sc[j]);
      std::size_t c = 1;
      while (!st.empty() && rank.back() == c) {
        cur = st.back() + cur; // earlier + later: the term order of the original left fold
        st.pop_back();
        rank.pop_back();
        c *= 2;
      }
      st.push_back(std::move(cur));
      rank.push_back(c);
    }

    if (st.empty()) return P(nsym);
    P acc = std::move(st.front()); // O(log n) leftovers, largest first
    for (std::size_t i = 1; i < st.size(); ++i)
      acc = acc + st[i];
    return acc;
  }

  /// @brief PHASE B, driver — fold every net, in parallel over the nets.
  ///
  /// Still net-parallel (unlike phase A), but that is fine here: the contraction is done, so a net's
  /// fold is proportional to its term count rather than to thousands of matrix products, and the
  /// merge in Codegen.m has already collapsed each net's repeated traces into one scalar apiece.
  template <class P, class TraceFn>
  std::vector<P> fold_nets(int nsym, const std::vector<std::vector<int>> &sidx,
                           const std::vector<std::vector<Cx>> &sc, const std::vector<P> &T, long nCache,
                           unsigned W, TraceFn &&trace)
  {
    std::vector<P> mp(sidx.size(), P(nsym));
    parallel_flat(static_cast<long>(sidx.size()), W, [&](long i) {
      const auto u = static_cast<std::size_t>(i);
      mp[u] = fold_net<P>(nsym, sidx[u], sc[u], T, nCache, trace);
    });
    return mp;
  }

} // namespace numtracer::numeric
