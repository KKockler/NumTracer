#include "../kernel.hh"

#include "../Zc.hh"

void Zc_integrator::get(double& dest, const double& p, const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA)
{
  const auto _nth = Zc_kernel<Regulator>::ntHoisted(k, ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA);
  integrator.get(dest, p, k, ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA, _nth[0], _nth[1], _nth[2], _nth[3], _nth[4], _nth[5]);
}