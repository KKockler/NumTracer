#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Gluon_condensate_plain_kernels.hh"

namespace DiFfRG
{
  class Gluon_condensate_plain_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflat, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflatDot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::gluon_condensate_plain::nenv) > 0 ? (DiFfRG::gluon_condensate_plain::nenv) : 1];
      DiFfRG::gluon_condensate_plain::fill(fenv, l1, cos1, p);
      const auto _interp1 = ZAflat(sqrt(powr<2>(l1)));
      const auto _interp2 = ZAflat(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = ZAflatDot(sqrt(powr<2>(l1)));
      return fma(-1., DiFfRG::gluon_condensate_plain::tr0(fenv) * powr<2>(_interp1) * _interp2 * _interp3, fma(DiFfRG::gluon_condensate_plain::tr1(fenv), powr<2>(_interp1) * _interp2 * _interp3, 0.));
    }

    static inline auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflat, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAflatDot)
    {
      return 0.;
    }
  };
}
using DiFfRG::Gluon_condensate_plain_kernel;