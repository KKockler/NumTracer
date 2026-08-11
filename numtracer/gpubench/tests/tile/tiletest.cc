// Validates DiFfRG::compute_divisible_tile against the cases tilewaste.py derives independently
// from the Kokkos sources. Host-only; no GPU needed.
#include <DiFfRG/common/kokkos.hh>
#include <cstdio>

using namespace DiFfRG;

static int failures = 0;

template <int dim>
void check(const char *label, device::array<size_t, dim> extents, device::array<size_t, dim> kokkos_tile,
           size_t budget, device::array<size_t, dim> want)
{
  const auto got = compute_divisible_tile<dim>(extents, kokkos_tile, budget);
  bool ok = true;
  size_t prod = 1, launched = 1, useful = 1;
  for (int i = 0; i < dim; ++i) {
    ok &= (got[i] == want[i]);
    prod *= got[i];
    launched *= ((extents[i] + got[i] - 1) / got[i]) * got[i];
    useful *= extents[i];
  }
  if (!ok) ++failures;
  printf("  [%s] %-28s got {", ok ? "ok " : "FAIL", label);
  for (int i = 0; i < dim; ++i)
    printf("%zu%s", got[i], i + 1 < dim ? "," : "");
  printf("} block=%zu waste=%.4fx", prod, double(launched) / double(useful));
  if (!ok) {
    printf("  WANT {");
    for (int i = 0; i < dim; ++i)
      printf("%zu%s", want[i], i + 1 < dim ? "," : "");
    printf("}");
  }
  printf("\n");
}

int main()
{
  printf("compute_divisible_tile:\n");

  // nf2 / with_mesons / with_mesons_3D rank-5 flows: the broken case. Kokkos picks tile 4 against
  // extent 6 at LB=128; the fix is {32,2,2,1,1}: still 128 threads (4 whole warps), zero masked lanes.
  check<5>("nf2 rank5 @128", {64, 32, 6, 6, 6}, {16, 2, 4, 1, 1}, 128, {32, 2, 2, 1, 1});

  // YangMills: orders are 8, the power-of-two tile already divides -> must be returned untouched.
  check<5>("ym rank5 @128 (clean)", {64, 32, 8, 8, 8}, {16, 2, 4, 1, 1}, 128, {16, 2, 4, 1, 1});

  // nf2 rank-4 flows are already clean at LB=128 -> untouched.
  check<4>("nf2 rank4 @128 (clean)", {64, 32, 6, 6}, {16, 2, 2, 2}, 128, {16, 2, 2, 2});

  // At LB=96 Kokkos already collapses to a divisible tile (this is why LB=96 looked like a win).
  check<5>("nf2 rank5 @96 (clean)", {64, 32, 6, 6, 6}, {16, 2, 2, 1, 1}, 64, {16, 2, 2, 1, 1});

  // No launch bounds: budget 256, Kokkos {16,2,4,2,1} still masks dim2.
  check<5>("nf2 rank5 @256", {64, 32, 6, 6, 6}, {16, 2, 4, 2, 1}, 256, {32, 2, 2, 2, 1});

  // finite-T: last axis is Matsubara (runtime length); waste is still on dim2.
  check<5>("fT rank5 @128, nmats=16", {64, 24, 6, 6, 16}, {16, 2, 4, 1, 1}, 128, {32, 2, 2, 1, 1});

  // Degenerate: a prime extent smaller than the tile must still divide (tile collapses to 1).
  check<3>("prime extent", {64, 32, 7}, {32, 2, 4}, 256, {64, 4, 1});

  printf("\n%s\n", failures == 0 ? "PASS" : "FAILURES");
  return failures != 0;
}
