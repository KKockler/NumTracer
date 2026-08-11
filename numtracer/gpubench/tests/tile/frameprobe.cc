// Stack-frame probe for DIFFRG_REF_ARGS.
//
// Reproduces the production QCD_Nf2 launch shape -- Integrator_p2_4D_3ang over {64,32,6,6,6} with
// *13* SplineInterpolator1D arguments, exactly what the generated kernels take -- so that
// `ptxas -v` reports the same per-thread stack frame the real flows carry.
//
// A100_NOTES.md 4.1 records a per-binary frame floor of 4000 B (Nf2, 13 interpolators) and 2416 B
// (YangMills, 7). sizeof(SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>)
// is 272 B, so the by-value tuple_cat in quadrature_integrator.hh accounts for 13*272 = 3536 B of
// that -- 89%. Compile this twice, -DDIFFRG_REF_ARGS=0 and =1, and compare the "stack frame" line.
//
// Build with ./run.sh --frame (or by hand with -Xptxas -v -cubin).
#include <DiFfRG/common/init.hh>
#include <DiFfRG/common/math.hh>
#include <DiFfRG/physics/integration/vacuum/integrator_p2_4d_3ang.hh>
#include <DiFfRG/physics/interpolation.hh>

#include <cstdio>
#include <vector>

using namespace DiFfRG;

using Interp = SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>;

// 13 interpolators, matching the nf2 generated kernel signature (ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4,
// ZAqbq7, dtZc, Zc, dtZA, ZA, zq, ZAAqbq1, ZAAqbq2). Each is read at a distinct argument so none
// can be optimised away.
class ProbeKernel13
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
  {
    constexpr size_t p_grid_size = 64;
    const std::array<size_t, 4> grid_size{32, 6, 6, 6};

    QuadratureProvider quadrature_provider;
    Integrator_p2_4D_3ang<4, double, ProbeKernel13, GPU_exec> integrator(quadrature_provider, grid_size, 100.);

    LogarithmicCoordinates1D<double> icoords(64, 1e-3, 1e2, 9.);
    std::vector<double> data(64);
    for (size_t i = 0; i < data.size(); ++i)
      data[i] = 1. + 0.01 * double(i);

    Interp interp(icoords);
    interp.update(data.data());

    LogarithmicCoordinates1D<float> coordinates(p_grid_size, 1e-3f, 1e2f, 9.f);
    std::vector<double> result(p_grid_size, 0.);
    integrator.map(result.data(), coordinates, interp, interp, interp, interp, interp, interp, interp, interp, interp,
                   interp, interp, interp, interp);
    Kokkos::fence();

    for (size_t i = 0; i < p_grid_size; ++i)
      printf("%3zu %.17e\n", i, result[i]);
  }
  return 0;
}
