// Timing harness for DIFFRG_DIVISIBLE_TILE x DIFFRG_REF_ARGS on whatever GPU is present.
//
// Same launch shape as frameprobe.cc (Integrator_p2_4D_3ang over {64,32,6,6,6} with 13 real
// interpolators), repeated so the measurement is not dominated by launch latency.
//
// The 4070 Laptop idles at 210 MHz and needs order a second of sustained load to reach boost, and
// that ramp has already faked a "1.3-1.5x" result in this project
// (docs/NUMTRACER_GPU_INVESTIGATION.md 5). So: a fixed-DURATION burn-in first, then the timed
// window, and bench.sh alternates the variants and reports ratios rather than absolutes.
#include <DiFfRG/common/init.hh>
#include <DiFfRG/common/math.hh>
#include <DiFfRG/physics/integration/vacuum/integrator_p2_4d_3ang.hh>
#include <DiFfRG/physics/interpolation.hh>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace DiFfRG;
using Interp = SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>;

class BenchKernel13
{
public:
  static KOKKOS_FORCEINLINE_FUNCTION double
  kernel(const double q, const double cos1, const double cos2, const double phi, const double p, const Interp &i0,
         const Interp &i1, const Interp &i2, const Interp &i3, const Interp &i4, const Interp &i5, const Interp &i6,
         const Interp &i7, const Interp &i8, const Interp &i9, const Interp &i10, const Interp &i11, const Interp &i12)
  {
    const double a = q * (1. + 0.5 * cos1) * (1. + 0.25 * cos2) * (1. + 0.125 * Kokkos::sin(phi));
    double acc = 1. / (1. + q * q) * (1. + 0.0625 * p);
    acc *= i0(a) + i1(a * 1.1) + i2(a * 1.2) + i3(a * 1.3) + i4(a * 1.4) + i5(a * 1.5) + i6(a * 1.6);
    acc *= i7(a * 1.7) + i8(a * 1.8) + i9(a * 1.9) + i10(a * 2.0) + i11(a * 2.1) + i12(a * 2.2);
    return acc;
  }
  static KOKKOS_FORCEINLINE_FUNCTION double constant(const double, const Interp &, const Interp &, const Interp &,
                                                     const Interp &, const Interp &, const Interp &, const Interp &,
                                                     const Interp &, const Interp &, const Interp &, const Interp &,
                                                     const Interp &, const Interp &)
  {
    return 0.;
  }
};

int main(int argc, char **argv)
{
  DiFfRG::Init(argc, argv);
  int rc = 0;
  {
    const int iters = (argc > 1) ? std::atoi(argv[1]) : 200;
    const double burn_seconds = (argc > 2) ? std::atof(argv[2]) : 3.0;

    constexpr size_t p_grid_size = 64;
    const std::array<size_t, 4> grid_size{32, 6, 6, 6};

    QuadratureProvider quadrature_provider;
    Integrator_p2_4D_3ang<4, double, BenchKernel13, GPU_exec> integrator(quadrature_provider, grid_size, 100.);

    LogarithmicCoordinates1D<double> icoords(64, 1e-3, 1e2, 9.);
    std::vector<double> data(64);
    for (size_t i = 0; i < data.size(); ++i)
      data[i] = 1. + 0.01 * double(i);
    Interp interp(icoords);
    interp.update(data.data());

    LogarithmicCoordinates1D<float> coordinates(p_grid_size, 1e-3f, 1e2f, 9.f);
    std::vector<double> result(p_grid_size, 0.);

    const auto once = [&]() {
      integrator.map(result.data(), coordinates, interp, interp, interp, interp, interp, interp, interp, interp,
                     interp, interp, interp, interp, interp);
      Kokkos::fence();
    };

    // Burn for a fixed wall duration so the clock is at boost before anything is timed.
    using clock = std::chrono::steady_clock;
    const auto burn_end = clock::now() + std::chrono::duration<double>(burn_seconds);
    long burned = 0;
    while (clock::now() < burn_end) {
      once();
      ++burned;
    }

    const auto t0 = clock::now();
    for (int i = 0; i < iters; ++i)
      once();
    const double elapsed = std::chrono::duration<double>(clock::now() - t0).count();

    // A checksum, so a variant that silently computes something else cannot look fast.
    double checksum = 0.;
    for (const double v : result)
      checksum += v;

    printf("iters=%d burn=%ld wall=%.6f per_iter_ms=%.6f checksum=%.17e\n", iters, burned, elapsed,
           1e3 * elapsed / iters, checksum);
  }
  return rc;
}
