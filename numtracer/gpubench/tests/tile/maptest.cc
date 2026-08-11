// A/B correctness check for the divisibility-aware MDRange tile.
//
// Runs the real Integrator_p2_4D_3ang::map() over the production QCD_Nf2 extents
// {p_grid=64, x_order=32, cos1=6, cos2=6, phi=6} -- the shape whose rank-5 tile Kokkos gets wrong
// -- and dumps every result. Compile this TU twice, once with -DDIFFRG_DIVISIBLE_TILE=0 and once
// with =1, and diff the dumps: re-tiling a parallel_for cannot change any arithmetic, so the two
// must be bit-identical. If they are not, the re-tiling is dropping or double-counting points.
//
// The kernel is deliberately asymmetric in every argument, so a lost or duplicated quadrature
// point cannot cancel out of the sum.
#include <DiFfRG/common/init.hh>
#include <DiFfRG/common/math.hh>
#include <DiFfRG/physics/integration/vacuum/integrator_p2_4d_3ang.hh>
#include <DiFfRG/physics/interpolation.hh>

#include <cstdio>
#include <vector>

using namespace DiFfRG;

class ProbeKernel
{
public:
  static KOKKOS_FORCEINLINE_FUNCTION double kernel(const double q, const double cos1, const double cos2,
                                                   const double phi, const double p)
  {
    // Asymmetric in each variable and bounded, so every point contributes distinguishably.
    return 1. / (1. + q * q) * (1. + 0.5 * cos1) * (1. + 0.25 * cos2 * cos2) * (1. + 0.125 * Kokkos::sin(phi)) *
           (1. + 0.0625 * p);
  }
  static KOKKOS_FORCEINLINE_FUNCTION double constant(const double) { return 0.; }
};

int main(int argc, char **argv)
{
  DiFfRG::Init(argc, argv);
  {
    // Production QCD_Nf2 shape: p_grid_size = 64 external points, {32, 6, 6, 6} quadrature orders.
    constexpr size_t p_grid_size = 64;
    const std::array<size_t, 4> grid_size{32, 6, 6, 6};

    QuadratureProvider quadrature_provider;
    Integrator_p2_4D_3ang<4, double, ProbeKernel, GPU_exec> integrator(quadrature_provider, grid_size, 100.);

    LogarithmicCoordinates1D<float> coordinates(p_grid_size, 1e-3f, 1e2f, 9.f); // bias as in nf2 parameter.json

    // Control: show the tile this build actually launches with, so a bit-identical result cannot be
    // mistaken for "the new code path never fired". Same construction as map() uses internally.
    {
      device::array<size_t, 5> start{0, 0, 0, 0, 0};
      device::array<size_t, 5> end{p_grid_size, grid_size[0], grid_size[1], grid_size[2], grid_size[3]};
      auto space = GPU_exec{};
      auto policy = make_kokkos_nd_range_divisible<5, GPU_exec>(space, start, end);
      size_t block = 1, launched = 1, useful = 1;
      fprintf(stderr, "DIFFRG_DIVISIBLE_TILE=%d  tile {", DIFFRG_DIVISIBLE_TILE);
      for (int i = 0; i < 5; ++i) {
        const size_t t = static_cast<size_t>(policy.m_tile[i]);
        fprintf(stderr, "%zu%s", t, i < 4 ? "," : "");
        block *= t;
        launched *= ((end[i] + t - 1) / t) * t;
        useful *= end[i];
      }
      fprintf(stderr, "} block=%zu waste=%.4fx\n", block, double(launched) / double(useful));
    }

    std::vector<double> result(p_grid_size, 0.);
    integrator.map(result.data(), coordinates);
    Kokkos::fence();

    for (size_t i = 0; i < p_grid_size; ++i)
      printf("%3zu %.17e\n", i, result[i]);
  }
  return 0; // DiFfRG::Init registers Kokkos finalization; calling it here double-finalizes
}
