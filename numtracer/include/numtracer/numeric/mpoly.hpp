/// @file mpoly.hpp
/// @brief Multivariate polynomial over an arbitrary user-symbol set, carrying surviving inverse
///        (`1/k²`) atoms in the monomial key. The arithmetic core of the **numeric contraction
///        backend**: contract a diagram as 4×4 spinor matrix products whose entries are
///        these polynomials, so γ/metric/projector stay numeric and only the user's symbolic
///        momentum data enters the result — bounded by the matrix structure, not the `(2n−1)!!`
///        sp-monomial blowup of the symbolic `reduce_product` path.
///
/// A `MPoly` is a map from monomial → complex coefficient. A monomial is
///   (exponent vector over the `nsym` user symbols, a sorted multiset of inverse-atom ids).
/// The atom ids ride along so a transverse projector's `INV(k) = 1/k²` factor is tracked exactly:
/// when a numerator monomial is divisible by an atom whose **denominator is a single monomial**
/// (e.g. a unit-direction loop `k²=l1²`), @ref divThroughMonomialAtoms cancels it; an atom whose
/// denominator is a genuine polynomial (a shifted line `k=l−q`, `k²` non-monomial) survives and is
/// later lowered to an `inv` env slot. This is the general bare-loop cancellation, not a special case.
///
/// Map-based (sizes are tens of monomials with frame inputs — the regime where this backend wins);
/// the symbol *meaning* lives only in the kernel's `fill`, so the engine is frame-agnostic.
///
/// This header is the polynomial type alone; the 4×4 spinor matrix of `MPoly` (`Mat4`) and the
/// γ/slash builders that consume it live in `numeric/spinor_mat.hpp`.
#pragma once

#include "numtracer/core/cx.hpp"
#include "numtracer/third_party/gch/small_vector.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

namespace numtracer::numeric
{

  using numtracer::Cx;

  /// Inline capacities for @ref Mono. `e` is ALWAYS exactly `nsym` long and `nsym` is the flow's env
  /// size (single digits to low tens), while `atoms` is usually empty and rarely more than a few — so
  /// with these buffers a monomial holds its whole payload inline and never touches the heap.
  ///
  /// This matters because `operator*` materialises all n·m product monomials: with `std::vector`
  /// members that was one heap allocation per product term, and profiling the generator run showed
  /// malloc/free at ~17% of total time (plus the pointer-chasing it forces on every sort comparison).
  /// Larger flows overflow these buffers gracefully — small_vector just falls back to the heap, i.e.
  /// exactly the old behaviour, so correctness never depends on the capacity being big enough.
  ///
  /// Exponents and atom ids are `int16_t`, not `int`. Both are small by construction — an exponent is
  /// a polynomial degree (single digits in practice) and an atom id indexes the flow's atom table — so
  /// 32 bits each is pure width. It matters because `std::sort` in @ref MPoly::from_scratch physically
  /// shuffles `pair<Mono, Cx>` objects and the monomial is the bulk of one: halving the exponent width
  /// halves what every swap moves and doubles how many monomials fit in cache.
  using MonoExpT = std::int16_t;
  using MonoAtomT = std::int16_t;

  inline constexpr unsigned kMonoExpInline = 24; ///< covers nsym for the flows in practice
  inline constexpr unsigned kMonoAtomInline = 8;

  using MonoExp = gch::small_vector<MonoExpT, kMonoExpInline>;
  using MonoAtoms = gch::small_vector<MonoAtomT, kMonoAtomInline>;

  struct Mono;
  /// Scratch list of (monomial, coeff) handed to @ref MPoly::from_scratch. Measured on a real flow,
  /// from_scratch is called tens of millions of times with a MEAN of ~4 terms — it is not one big
  /// sort, it is a flood of tiny ones — so a heap-backed scratch vector paid an allocation per call.
  /// Inline storage covers the common case; larger products fall back to the heap as before.
  inline constexpr unsigned kMPolyScratchInline = 8;

  /// @brief A monomial: exponents over the user symbols plus a sorted multiset of inverse-atom ids.
  struct Mono {
    MonoExp e;       ///< length nsym, exponent of each user symbol
    MonoAtoms atoms; ///< sorted (with multiplicity) ids of surviving `1/k²` atoms
    bool operator<(const Mono &o) const
    {
      if (e != o.e) return e < o.e;
      return atoms < o.atoms;
    }
    bool operator==(const Mono &o) const { return e == o.e && atoms == o.atoms; }
  };

  using MPolyScratch = gch::small_vector<std::pair<Mono, Cx>, kMPolyScratchInline>;

  // The `nsym`-carrying construction API is closed behind these friends (see @ref LorentzEnv): every
  // polynomial in one trace must share an `nsym`, so the factories that bake it in are private and the
  // ONLY sanctioned way to mint one is a @ref LorentzEnv (which holds a fixed `nsym`). @ref MPolyFactory
  // is a tiny internal attorney that re-exposes the same factories to the trusted cross-header engine
  // code (contraction / trace-fold), which already threads a single `nsym` and must not route through a
  // user-facing env. Only names are needed here; the definitions live below / in `numeric/env.hpp`.
  class LorentzEnv;
  struct MPolyFactory;

  /// @brief Multivariate polynomial: monomial → complex coefficient, over `nsym` user symbols.
  ///
  /// Storage is a **sorted, like-terms-combined `std::vector`** (not a `std::map`): the per-insert
  /// red-black-tree overhead of a map dominates for the large (10⁴+ monomial) pure-gauge polynomials,
  /// so `operator*` collects all `n·m` products into a scratch vector and **sort-collects once**
  /// (O(nm log nm), no per-term tree churn) and `operator+` is a linear merge. `t` stays sorted by
  /// @ref Mono and carries no zero coefficients, so iteration order is deterministic (reproducible
  /// kernel) and equality/lookup are binary searches.
  struct MPoly {
    int nsym = 0;
    /// Heap-backed on purpose. Giving `t` inline storage (capacity 4) was measured: it bought ~2% of
    /// run time but pushed sizeof(Mat4) from ~0.5 KB to ~8.7 KB and peak RSS up ~20% — a bad trade on
    /// the dense flows, which are memory-bound before they are time-bound.
    std::vector<std::pair<Mono, Cx>> t; ///< sorted by Mono, like terms combined, no zeros

    // Sanctioned construction paths (see the note above @ref MPoly). The env and the internal attorney
    // reach the nsym-taking ctor/factories; the in-header arithmetic operators keep constructing result
    // polynomials directly (they already have a definite `nsym` from their operands), so their bodies
    // stay byte-for-byte unchanged — the hot path is untouched.
    friend class LorentzEnv;
    friend struct MPolyFactory;
    friend MPoly operator+(const MPoly &a, const MPoly &b);
    friend MPoly operator-(const MPoly &a, const MPoly &b);
    friend MPoly operator*(const MPoly &a, const MPoly &b);
    friend MPoly divThroughMonomialAtoms(const MPoly &p, const std::vector<MPoly> &atomDen);
    friend MPoly reduce_units(const MPoly &p, const std::vector<std::vector<int>> &groups);

    MPoly() = default;

  private:
    // Bare-`nsym` construction — reachable only through @ref LorentzEnv / @ref MPolyFactory (friends).
    // Making these private is what turns "every MPoly in one trace shares an nsym" from a debug-time
    // assert into a compile-time guarantee: outside the sanctioned env you cannot mint a non-empty MPoly
    // with a hand-picked nsym. The empty default ctor above stays public (an nsym==0 zero used by the
    // operator short-circuits and by std::array/std::vector default members).
    explicit MPoly(int ns) : nsym(ns) {}

    /// Build from an unsorted scratch list of (monomial, coeff): sort then combine adjacent equals.
    static MPoly from_scratch(int ns, MPolyScratch s)
    {
      MPoly p(ns);
      std::sort(s.begin(), s.end(), [](const auto &a, const auto &b) { return a.first < b.first; });
      p.t.reserve(s.size());
      for (auto &kv : s) {
        if (kv.second.re == 0 && kv.second.im == 0) continue;
        if (!p.t.empty() && p.t.back().first == kv.first) {
          p.t.back().second = p.t.back().second + kv.second;
          if (p.t.back().second.re == 0 && p.t.back().second.im == 0) p.t.pop_back();
        } else {
          p.t.push_back(std::move(kv));
        }
      }
      return p;
    }

    static MPoly constant(int ns, Cx c)
    {
      MPoly p(ns);
      if (!(c.re == 0 && c.im == 0)) p.t.push_back({Mono{MonoExp(ns, 0), {}}, c});
      return p;
    }
    /// The i-th user symbol (coefficient 1).
    static MPoly var(int ns, int i)
    {
      MPoly p(ns);
      MonoExp e(ns, 0);
      e[i] = 1;
      p.t.push_back({Mono{std::move(e), {}}, Cx{1, 0}});
      return p;
    }
    /// A single monomial `c · ∏ x_k^{e_k}` (no inverse atoms). Used by the generated component table.
    /// Takes a `std::vector` because the GENERATED component table hands one over as a braced list;
    /// it is converted into the monomial's inline storage here (once per table entry, not on any hot
    /// path).
    static MPoly mono(int ns, const std::vector<int> &e, Cx c)
    {
      MPoly p(ns);
      if (!(c.re == 0 && c.im == 0)) p.t.push_back({Mono{MonoExp(e.begin(), e.end()), {}}, c});
      return p;
    }
    /// A bare inverse atom `1/D` (atom id `aid`), coefficient 1.
    static MPoly atom(int ns, int aid)
    {
      MPoly p(ns);
      p.t.push_back({Mono{MonoExp(ns, 0), MonoAtoms{static_cast<MonoAtomT>(aid)}}, Cx{1, 0}});
      return p;
    }

  public:
    int size() const { return t.size(); }
    bool empty() const { return t.empty(); }

    /// Insert/accumulate one term (keeps `t` sorted). O(n) shift — used only on the incremental paths
    /// (component builders); the hot `operator*`/`operator+`/reductions build whole vectors at once.
    void addTerm(const Mono &m, Cx c)
    {
      if (c.re == 0 && c.im == 0) return;
      auto it = std::lower_bound(t.begin(), t.end(), m,
                                 [](const std::pair<Mono, Cx> &a, const Mono &k) { return a.first < k; });
      if (it != t.end() && it->first == m) {
        it->second = it->second + c;
        if (it->second.re == 0 && it->second.im == 0) t.erase(it);
      } else {
        t.insert(it, {m, c});
      }
    }
  };

  /// @brief Internal attorney re-exposing the private @ref MPoly factories to the trusted cross-header
  ///        engine code (spinor matrices, contraction, trace-fold). Those functions already carry one
  ///        definite `nsym` and must not depend on the user-facing @ref LorentzEnv, but they cannot be
  ///        friended by name (templates / DiracNet-heavy signatures across headers), so the friend
  ///        surface is localised to this one struct. NOT part of the public API — call sites outside the
  ///        engine construct polynomials through @ref LorentzEnv.
  struct MPolyFactory {
    static MPoly zero(int ns) { return MPoly(ns); }
    static MPoly constant(int ns, Cx c) { return MPoly::constant(ns, c); }
    static MPoly atom(int ns, int aid) { return MPoly::atom(ns, aid); }
    static MPoly from_scratch(int ns, MPolyScratch s) { return MPoly::from_scratch(ns, std::move(s)); }
  };

  inline MPoly operator+(const MPoly &a, const MPoly &b)
  {
    if (a.t.empty()) return b;
    if (b.t.empty()) return a;
    // Both operands carry terms, so their symbol spaces must agree: the merge below walks the two
    // exponent vectors slot-for-slot, and a mismatch would read out of bounds (a default-constructed
    // zero with nsym == 0 is allowed — it is caught by the empty checks above). Debug-only; compiles
    // out under NDEBUG.
    assert(a.nsym == b.nsym);
    const int ns = a.nsym ? a.nsym : b.nsym;
    MPoly r(ns);
    r.t.reserve(a.t.size() + b.t.size());
    std::size_t i = 0, j = 0;
    while (i < a.t.size() && j < b.t.size()) {
      if (a.t[i].first < b.t[j].first)
        r.t.push_back(a.t[i++]);
      else if (b.t[j].first < a.t[i].first)
        r.t.push_back(b.t[j++]);
      else {
        Cx s = a.t[i].second + b.t[j].second;
        if (!(s.re == 0 && s.im == 0)) r.t.push_back({a.t[i].first, s});
        ++i;
        ++j;
      }
    }
    while (i < a.t.size())
      r.t.push_back(a.t[i++]);
    while (j < b.t.size())
      r.t.push_back(b.t[j++]);
    return r;
  }
  inline MPoly operator-(const MPoly &a, const MPoly &b)
  {
    assert(a.t.empty() || b.t.empty() || a.nsym == b.nsym); // see operator+; debug-only
    MPoly nb(b.nsym);
    nb.t.reserve(b.t.size());
    for (const auto &kv : b.t)
      nb.t.push_back({kv.first, Cx{-kv.second.re, -kv.second.im}});
    return a + nb;
  }
  inline MPoly operator*(const MPoly &a, const MPoly &b)
  {
    const int ns = a.nsym ? a.nsym : b.nsym;
    if (a.t.empty() || b.t.empty()) return MPoly(ns);
    assert(a.nsym == b.nsym); // both carry terms ⇒ symbol spaces must match; see operator+ (debug-only)
    MPolyScratch s;
    s.reserve(a.t.size() * b.t.size());
    // Build each product monomial IN PLACE. The old code kept a scratch `e` and copied it into the
    // Mono, which (with a heap-backed exponent vector) was an allocation per product term.
    for (const auto &[ma, ca] : a.t)
      for (const auto &[mb, cb] : b.t) {
        auto &slot = s.emplace_back(Mono{}, ca * cb);
        Mono &m = slot.first;
        m.e.resize(static_cast<std::size_t>(ns));
        for (int k = 0; k < ns; ++k)
          m.e[k] = ma.e[k] + mb.e[k];
        // Only ask for capacity when the merge would actually overflow the inline buffer. Most
        // monomials carry NO atoms at all, and an unconditional reserve() here cost ~7% of the run.
        const std::size_t nat = ma.atoms.size() + mb.atoms.size();
        if (nat > kMonoAtomInline) m.atoms.reserve(nat);
        std::size_t i = 0, j = 0; // merge the two sorted atom multisets
        while (i < ma.atoms.size() && j < mb.atoms.size())
          m.atoms.push_back(ma.atoms[i] <= mb.atoms[j] ? ma.atoms[i++] : mb.atoms[j++]);
        while (i < ma.atoms.size())
          m.atoms.push_back(ma.atoms[i++]);
        while (j < mb.atoms.size())
          m.atoms.push_back(mb.atoms[j++]);
      }
    return MPoly::from_scratch(ns, std::move(s));
  }

  /// @brief Cancel each numerator monomial against any atom whose denominator is a single monomial.
  ///
  /// `atomDen[aid]` is the atom's denominator polynomial `k²`. If it is a single monomial
  /// `c·∏ x^d`, then one power of the atom `1/D` cancels whenever the numerator exponents dominate
  /// `d` componentwise: subtract `d`, divide the coefficient by `c`, drop one atom instance. Repeat
  /// per atom occurrence. Atoms with a non-monomial denominator (shifted lines) are left intact and
  /// survive into the lowering as `inv` env slots. Value-preserving and frame-agnostic.
  inline MPoly divThroughMonomialAtoms(const MPoly &p, const std::vector<MPoly> &atomDen)
  {
    MPolyScratch out;
    out.reserve(p.t.size());
    for (const auto &[m, c] : p.t) {
      MonoExp e = m.e;
      Cx coeff = c;
      MonoAtoms keep;
      keep.reserve(m.atoms.size());
      for (const auto &aid : m.atoms) {
        bool cancelled = false;
        if (aid >= 0 && aid < (int)atomDen.size()) {
          const MPoly &D = atomDen[aid];
          if (D.t.size() == 1) { // monomial denominator → may cancel
            const auto &dm = *D.t.begin();
            const MonoExp &d = dm.first.e;
            const Cx dc = dm.second;
            bool dominates = dm.first.atoms.empty();
            for (std::size_t k = 0; dominates && k < e.size(); ++k)
              if (e[k] < d[k]) dominates = false;
            if (dominates) {
              for (std::size_t k = 0; k < e.size(); ++k)
                e[k] -= d[k];
              // coeff /= dc   (complex division: z / w = z·conj(w) / |w|²)
              const double den = dc.re * dc.re + dc.im * dc.im;
              assert(den != 0.0); // a stored monomial denominator never has a zero coefficient
              coeff = Cx{(coeff.re * dc.re + coeff.im * dc.im) / den, (coeff.im * dc.re - coeff.re * dc.im) / den};
              cancelled = true;
            }
          }
        }
        if (!cancelled) keep.push_back(aid);
      }
      // keep stays sorted because m.atoms was sorted and we only dropped elements
      if (!(coeff.re == 0 && coeff.im == 0)) out.push_back({Mono{std::move(e), std::move(keep)}, coeff});
    }
    return MPoly::from_scratch(p.nsym, std::move(out));
  }

  /// @brief Apply unit-vector constraints `Σ_μ Uμ² = 1` to reduce each unit group's LAST component to
  ///        power ≤ 1: `U_last² → 1 − Σ_{μ<last} Uμ²`. A `group` is the list of symbol indices of one
  ///        unit direction's components (e.g. the loop direction `{U0,U1,U2}` with `U0²+U1²+U2²=1`, or a
  ///        `{cos,sin}` pair). This is the general form of `sin²→1−cos²` and is the crux of the numeric
  ///        backend's compactness: it (a) collapses the bare-loop denominator `Σ(l1·Uμ)² = l1²·ΣUμ²` to
  ///        the **monomial `l1²`** so @ref divThroughMonomialAtoms cancels it (like inv's `rel`), and
  ///        (b) collapses the `U·U` factors transverse projectors generate. Value-preserving;
  ///        terminates because each rewrite strictly lowers the last component's exponent.
  inline MPoly reduce_units(const MPoly &p, const std::vector<std::vector<int>> &groups)
  {
    if (groups.empty()) return p;
    // every group entry is a symbol index, so it must address a valid component slot `e[idx]`
    for (const auto &g : groups)
      for ([[maybe_unused]] int idx : g)
        assert(idx >= 0 && idx < p.nsym);
    MPolyScratch out;
    std::vector<std::tuple<MonoExp, MonoAtoms, Cx>> work;
    for (const auto &[m, c] : p.t)
      work.push_back({m.e, m.atoms, c});
    while (!work.empty()) {
      auto [e, atoms, c] = std::move(work.back());
      work.pop_back();
      // find a unit group whose LAST component still has power >= 2 (the rewrite target)
      int groupIdx = -1;
      for (int gIdx = 0; gIdx < (int)groups.size(); ++gIdx)
        if (!groups[gIdx].empty() && e[groups[gIdx].back()] >= 2) {
          groupIdx = gIdx;
          break;
        }
      if (groupIdx < 0) {
        out.push_back({Mono{std::move(e), std::move(atoms)}, c});
        continue;
      }
      const std::vector<int> &group = groups[groupIdx];
      const int last = group.back();
      MonoExp base = e;
      base[last] -= 2;                                     // U_last^2 -> 1 - Σ_{μ<last} Uμ^2
      work.push_back({base, atoms, c});                    // the "+1" branch
      for (std::size_t i = 0; i + 1 < group.size(); ++i) { // the "-Uμ^2" branches
        MonoExp shifted = base;
        shifted[group[i]] += 2;
        work.push_back({std::move(shifted), atoms, Cx{-c.re, -c.im}});
      }
    }
    return MPoly::from_scratch(p.nsym, std::move(out));
  }

  /// @brief Numeric evaluation (validation only). `x[i]` = user symbol i; `atomVal[aid]` = value of
  ///        `1/D_aid` (the caller supplies the reciprocal already evaluated).
  inline Cx eval(const MPoly &p, const std::vector<double> &x, const std::vector<double> &atomVal)
  {
    Cx s{0, 0};
    for (const auto &[m, c] : p.t) {
      double mon = 1.0;
      for (int k = 0; k < p.nsym; ++k)
        for (int j = 0; j < m.e[k]; ++j)
          mon *= x[k];
      for (int aid : m.atoms)
        mon *= atomVal[aid];
      s = s + Cx{c.re * mon, c.im * mon};
    }
    return s;
  }

} // namespace numtracer::numeric
