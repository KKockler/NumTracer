#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Zc_kernel
  {
    public:
    using Regulator = REG;

    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zc_ym::nenv) > 0 ? (DiFfRG::zc_ym::nenv) : 1];
      DiFfRG::zc_ym::fill(fenv, l1, cos1, p);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _ix1 = ZA.index(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp4 = ZA.at(_ix1);
      const auto _interp5 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _ix0 = Zc.index(l1);
      const auto _interp6 = Zc.at(_ix0);
      const auto _interp7 = ZA.at(_ix0);
      const auto _interp8 = Zc.at(_ix1);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp7 * powr<2>(l1) + _interp1 * nthk0);
      const auto _den3 = powr<-2>(_interp6 * powr<2>(l1) + _interp1 * nthk1);
      const auto _den4 = powr<-1>(_interp3 * nthk0 + _interp4 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp3 * nthk1 + _interp8 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      return fma(-0.125, _cse1 * _den3 * _den4 * DiFfRG::zc_ym::tr0(fenv) * powr<2>(_interp5) * (_interp2 * nthk1 + _interp1 * (nthk2 + 50. * (-nthk1 + nthk3))), fma(0.125, _cse1 * _den2 * _den5 * DiFfRG::zc_ym::tr1(fenv) * powr<2>(_interp5) * (_interp2 * nthk0 + _interp1 * (nthk4 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk5))), 0.));
    }

    static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    {
      return 0.;
    }

    static device::array<double, 6> ntHoisted(const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      return {{ZA.CPU()(pow(1. + powr<6>(k),0.16666666666666666667)),
          Zc.CPU()(k),
          dtZc.CPU()(k),
          Zc.CPU()(1.02 * k),
          dtZA.CPU()(pow(1. + powr<6>(k),0.16666666666666666667)),
          ZA.CPU()(1.02 * pow(1. + powr<6>(k),0.16666666666666666667))}};
    }
    private:
    static KOKKOS_INLINE_FUNCTION auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static KOKKOS_INLINE_FUNCTION auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static KOKKOS_INLINE_FUNCTION auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static KOKKOS_INLINE_FUNCTION auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static KOKKOS_INLINE_FUNCTION auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static KOKKOS_INLINE_FUNCTION auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using DiFfRG::Zc_kernel;