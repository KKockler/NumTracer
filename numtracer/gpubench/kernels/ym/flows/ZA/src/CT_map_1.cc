#include "../kernel.hh"

#include "../ZA.hh"

DiFfRG::GPU_exec ZA_integrator::map(double* dest, const FocusedLogCoordinates1D<double>& coordinates, const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA)
{
  const auto _nth = ZA_kernel<Regulator>::ntHoisted(k, ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA);
  return integrator.map(dest, coordinates, k, ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA, _nth[0], _nth[1], _nth[2], _nth[3], _nth[4], _nth[5]);
}