/// @file real_cse.hpp
/// @brief The runtime real-SSA CSE builder driven by `codegen/lower.hpp::horner`.
///
/// Lowering a contracted diagram polynomial to straight-line code needs a place to accumulate
/// real arithmetic instructions while sharing every repeated subexpression exactly once. That is
/// @ref numtracer::network::rdetail::RBuilder: a growable array of real SSA instructions
/// (@ref RInstr) with hash-consed value numbering — `find_or_add` returns the existing slot for
/// an identical instruction instead of appending a duplicate. The `rconst`/`rvar`/`rneg`/`rmul`/
/// `radd`/`rsub` helpers layer the cheap algebraic folds on top (`×0`/`×1`/`×(-1)`, `+0`,
/// `-(-x)`, commutative operand ordering), so the emitted program is already free of trivial ops
/// before it reaches the C++ printer in `codegen/gen.hpp`.
///
/// This is a **build-time / generation-time** structure: the generator contracts each diagram
/// numerically, Horner-factors the resulting monomials into an `RBuilder`, and prints the
/// straight-line result. It is not part of the consumer-side kernel.
#pragma once

#include <bit>
#include <cstdint>
#include <vector>

namespace numtracer::network
{

  /// @brief Opcode tag for a real-SSA instruction.
  enum ROp : int {
    RCONST, ///< Real constant (value in @ref RInstr::k).
    RVAR,   ///< Real variable (environment index in @ref RInstr::a).
    RADD,   ///< `a + b`.
    RSUB,   ///< `a - b`.
    RMUL,   ///< `a * b`.
    RNEG    ///< `-a`.
  };
  /// @brief One instruction of the lowered real SSA.
  struct RInstr {
    int op = RCONST; ///< Opcode (one of @ref ROp).
    int a = -1;      ///< First operand slot; for `RVAR` this is the environment id.
    int b = -1;      ///< Second operand slot (binary ops only).
    double k = 0;    ///< Constant value, for `RCONST`.
  };

  namespace rdetail
  {
    /// @brief splitmix64 finalizer — scrambles a 64-bit value (local copy; any hash serves, the
    ///        emitted SSA depends only on insertion order, not on the hash).
    constexpr std::uint64_t mix64(std::uint64_t x)
    {
      x ^= x >> 30;
      x *= 0xbf58476d1ce4e5b9ULL;
      x ^= x >> 27;
      x *= 0x94d049bb133111ebULL;
      x ^= x >> 31;
      return x;
    }
    /// @brief Combine two hashes into one (boost-style mix; order-sensitive).
    constexpr std::uint64_t hcomb(std::uint64_t a, std::uint64_t b)
    {
      return mix64(a ^ (b + 0x9e3779b97f4a7c15ULL + (a << 6) + (a >> 2)));
    }

    // ---- hash-indexed buffer for the real SSA: no fixed capacity, O(1) amortised dedup -----
    // A dynamically growing instruction array so the real lowering has no capacity bound *and*
    // value-numbering stays O(w) (not O(w²)) as the real SSA grows.
    struct RBuilder {
      std::vector<RInstr> ins; ///< Real instructions, in emission order.
      std::vector<int> bucket; ///< Open-addressed index: slot, or -1 if empty.
      std::size_t mask = 0;    ///< `bucket.size()-1` (power of two); 0 while empty.

      static constexpr std::uint64_t ihash(const RInstr &e)
      {
        const std::uint64_t h = hcomb(hcomb(e.op, e.a), e.b);
        return hcomb(h, std::bit_cast<std::uint64_t>(e.k));
      }
      static constexpr bool ieq(const RInstr &a, const RInstr &b)
      {
        return a.op == b.op && a.a == b.a && a.b == b.b && a.k == b.k;
      }
      constexpr void rehash(std::size_t cap)
      {
        bucket.assign(cap, -1);
        mask = cap - 1;
        for (int s = 0; s < (int)ins.size(); ++s) {
          std::size_t p = ihash(ins[s]) & mask;
          while (bucket[p] != -1)
            p = (p + 1) & mask;
          bucket[p] = s;
        }
      }
      /// @brief Append an instruction, or reuse an identical existing one (value numbering).
      constexpr int find_or_add(RInstr e)
      {
        if (mask) {
          std::size_t p = ihash(e) & mask;
          while (bucket[p] != -1) {
            if (ieq(ins[bucket[p]], e)) return bucket[p];
            p = (p + 1) & mask;
          }
        }
        if ((ins.size() + 1) * 10 >= (mask + 1) * 7) rehash(mask == 0 ? 16 : (mask + 1) * 2);
        const int s = static_cast<int>(ins.size());
        ins.push_back(e);
        std::size_t p = ihash(e) & mask;
        while (bucket[p] != -1)
          p = (p + 1) & mask;
        bucket[p] = s;
        return s;
      }
    };

    /// @brief Emit a real constant, mapping `0` to the structural-zero sentinel (`-1`).
    constexpr int rconst(RBuilder &w, double k) { return k == 0.0 ? -1 : w.find_or_add({RCONST, -1, -1, k}); }
    /// @brief Emit a real variable (environment lookup).
    constexpr int rvar(RBuilder &w, int id) { return w.find_or_add({RVAR, id, -1, 0}); }
    /// @brief Emit a negation, folding `-(const)` and cancelling `-(-x)`.
    constexpr int rneg(RBuilder &w, int y)
    {
      if (y < 0) return -1;
      if (w.ins[y].op == RCONST) return w.find_or_add({RCONST, -1, -1, -w.ins[y].k});
      if (w.ins[y].op == RNEG) return w.ins[y].a;
      return w.find_or_add({RNEG, y, -1, 0});
    }
    /// @brief Emit a multiply, folding `×0`, `×1`, `×(-1)` and canonicalising operand order.
    constexpr int rmul(RBuilder &w, int x, int y)
    {
      if (x < 0 || y < 0) return -1;
      if (w.ins[x].op == RCONST) {
        if (w.ins[x].k == 1.0) return y;
        if (w.ins[x].k == -1.0) return rneg(w, y);
      }
      if (w.ins[y].op == RCONST) {
        if (w.ins[y].k == 1.0) return x;
        if (w.ins[y].k == -1.0) return rneg(w, x);
      }
      return w.find_or_add({RMUL, x < y ? x : y, x < y ? y : x, 0});
    }
    /// @brief Emit an addition, dropping structural-zero operands and canonicalising order.
    constexpr int radd(RBuilder &w, int x, int y)
    {
      if (x < 0) return y;
      if (y < 0) return x;
      return w.find_or_add({RADD, x < y ? x : y, x < y ? y : x, 0});
    }
    /// @brief Emit a subtraction, handling structural-zero operands.
    constexpr int rsub(RBuilder &w, int x, int y)
    {
      if (y < 0) return x;
      if (x < 0) return rneg(w, y);
      return w.find_or_add({RSUB, x, y, 0});
    }
  } // namespace rdetail

} // namespace numtracer::network
