/// @file lower.hpp
/// @brief The Horner core that lowers a set of monomials (@ref LMono) to a real straight-line
///        program through the shared CSE back-end.
///
/// A contracted diagram polynomial is a set of monomials, each a real coefficient times a product
/// of `var(envId)^exp` factors (the env id names a fundamental symbol — a scalar product, an
/// inverse propagator, a dressing call, or a raw kernel argument). Evaluating that fast is the job
/// of hash-consed value-numbering (@ref numtracer::network::rdetail::RBuilder — the CSE) plus greedy
/// multivariate Horner factoring: @ref horner emits the monomial set into an `RBuilder`, sharing
/// powers / products / factors, and returns the result slot. The numeric backend
/// (@ref numtracer::numeric::to_genprog) builds the @ref LMono set and drives this via
/// @ref numtracer::network::gdetail::best_into in `codegen/gen.hpp`.
#pragma once

#include "numtracer/codegen/real_cse.hpp"

#include <algorithm>
#include <climits>
#include <utility>
#include <vector>

namespace numtracer::network {

/// @brief A monomial of the lowered polynomial: `c * prod(var(envId)^exp)`, `vp` sorted by envId.
struct LMono {
  double c = 0;
  std::vector<std::pair<int, int>> vp; ///< (envId, exponent)
};

/// @brief Pick the Horner pivot: the variable occurring in the most terms (a linear tally — no
///        std::map, to stay constexpr), and the lowest exponent `pivotExp` it appears with (the most
///        we can factor out). Returns `{pivot envId, pivotExp}`.
constexpr std::pair<int, int> choose_pivot(const std::vector<LMono> &terms) {
  std::vector<std::pair<int, int>> varCounts; // (envId, #terms containing it)
  for (const LMono &m : terms)
    for (auto [id, e] : m.vp) {
      bool found = false;
      for (auto &vc : varCounts) if (vc.first == id) { ++vc.second; found = true; break; }
      if (!found) varCounts.push_back({id, 1});
    }
  int pivot = -1, bestCount = -1;
  for (auto &vc : varCounts)
    if (vc.second > bestCount) { bestCount = vc.second; pivot = vc.first; }
  int pivotExp = INT_MAX;
  for (const LMono &m : terms)
    for (auto [id, ex] : m.vp)
      if (id == pivot) pivotExp = std::min(pivotExp, ex);
  return {pivot, pivotExp};
}

/// @brief Partition @p terms into `with` (those containing `pivot`, with `pivot^pivotExp` divided
///        out) and `without` (the rest). Returns `{with, without}`.
constexpr std::pair<std::vector<LMono>, std::vector<LMono>>
partition_pivot(std::vector<LMono> terms, int pivot, int pivotExp) {
  std::vector<LMono> with, without;
  for (LMono &m : terms) {
    bool has = false;
    for (auto [id, ex] : m.vp)
      if (id == pivot) { has = true; break; }
    if (has) {
      LMono reduced;
      reduced.c = m.c;
      for (auto [id, ex] : m.vp) {
        if (id == pivot) { if (ex > pivotExp) reduced.vp.push_back({id, ex - pivotExp}); }
        else reduced.vp.push_back({id, ex});
      }
      with.push_back(std::move(reduced));
    } else {
      without.push_back(std::move(m));
    }
  }
  return {std::move(with), std::move(without)};
}

/// @brief Greedy multivariate Horner of a monomial set, emitted through the real value-numbering
///        builder (so shared powers / products / factors are CSE'd). Returns the result slot.
constexpr int horner(rdetail::RBuilder &w, std::vector<LMono> terms) {
  using namespace rdetail;
  // ---- base case: all monomials are constant → emit their sum -----------------------------
  bool allconst = true;
  for (const LMono &m : terms)
    if (!m.vp.empty()) { allconst = false; break; }
  if (allconst) {
    double s = 0;
    for (const LMono &m : terms) s += m.c;
    return rconst(w, s);
  }
  const auto [pivot, pivotExp] = choose_pivot(terms);
  auto [with, without] = partition_pivot(std::move(terms), pivot, pivotExp);
  // ---- recurse on both parts, then combine: pivot^pivotExp * horner(with) + horner(without) -----
  const int withSlot = horner(w, std::move(with));
  const int withoutSlot = without.empty() ? -1 : horner(w, std::move(without));
  const int pivotVar = rvar(w, pivot);
  int pivotPow = pivotVar;
  for (int e = 1; e < pivotExp; ++e) pivotPow = rmul(w, pivotPow, pivotVar);
  const int factored = rmul(w, pivotPow, withSlot);
  return radd(w, factored, withoutSlot);
}

} // namespace numtracer::network
