#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Gluon_condensate_diag_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Gluon_condensate_diag_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA8, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAoff, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3dot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA8dot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAoffdot)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::gluon_condensate_diag::nenv) > 0 ? (DiFfRG::gluon_condensate_diag::nenv) : 1];
      DiFfRG::gluon_condensate_diag::fill(fenv, l1, cos1, p);
      const auto _interp1 = ZA3(sqrt(powr<2>(l1)));
      const auto _interp2 = ZA3(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = ZA3dot(sqrt(powr<2>(l1)));
      const auto _interp4 = ZA8(sqrt(powr<2>(l1)));
      const auto _interp5 = ZA8(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA8dot(sqrt(powr<2>(l1)));
      const auto _interp7 = ZAoff(sqrt(powr<2>(l1)));
      const auto _interp8 = ZAoff(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZAoffdot(sqrt(powr<2>(l1)));
      return 6. * fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp1) * _interp2 * _interp3, fma(DiFfRG::gluon_condensate_diag::tr0(fenv), powr<2>(_interp4) * _interp5 * _interp6, fma(6., DiFfRG::gluon_condensate_diag::tr0(fenv) * powr<2>(_interp7) * _interp8 * _interp9, 0.)));
    }

    static inline auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA8, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAoff, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3dot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA8dot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAoffdot)
    {
      return 0.;
    }private: static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using DiFfRG::Gluon_condensate_diag_kernel;