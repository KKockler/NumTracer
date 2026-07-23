#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Flavour_split_num_kernels.hh"

namespace DiFfRG
{
  class Flavour_split_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gu, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gd, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GuDot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GdDot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::flavour_split_num::nenv) > 0 ? (DiFfRG::flavour_split_num::nenv) : 1];
      DiFfRG::flavour_split_num::fill(fenv, l1, cos1, p);
      const auto _interp1 = Gd(sqrt(powr<2>(l1)));
      const auto _interp2 = Gd(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = GdDot(sqrt(powr<2>(l1)));
      const auto _interp4 = Gu(sqrt(powr<2>(l1)));
      const auto _interp5 = Gu(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = GuDot(sqrt(powr<2>(l1)));
      return -2. * fma(DiFfRG::flavour_split_num::tr1(fenv), powr<2>(_interp1) * _interp2 * _interp3, fma(DiFfRG::flavour_split_num::tr0(fenv), powr<2>(_interp4) * _interp5 * _interp6, 0.));
    }

    static inline auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gu, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gd, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GuDot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GdDot)
    {
      return 0.;
    }
  };
}
using DiFfRG::Flavour_split_num_kernel;