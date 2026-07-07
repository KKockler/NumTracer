#pragma once

#include "shim.hpp"

namespace DiFfRG {
  template<typename _Regulator>
  class Flavour_split_form_kernel
  {
    public:
    using Regulator = _Regulator;

    static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gu, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gd, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GuDot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GdDot)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;const auto _interp1 = Gd(l1);
      const auto _interp2 = Gd(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp3 = GdDot(l1);
      const auto _interp4 = Gu(l1);
      const auto _interp5 = Gu(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = GuDot(l1);
      const auto _cse1 = -cos1 * p;
      const auto _cse2 = _cse1 + l1;
      return -24. * fma(_cse2, powr<2>(_interp1) * _interp2 * _interp3 * powr<3>(l1), fma(_cse2, powr<2>(_interp4) * _interp5 * _interp6 * powr<3>(l1), 0.));
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gu, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Gd, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GuDot, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& GdDot)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;
      return 0.;
    }private: static KOKKOS_FORCEINLINE_FUNCTION auto RB(const auto& k2, const auto& p2)
    {
      return Regulator::RB(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto RF(const auto& k2, const auto& p2)
    {
      return Regulator::RF(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto RBdot(const auto& k2, const auto& p2)
    {
      return Regulator::RBdot(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto RFdot(const auto& k2, const auto& p2)
    {
      return Regulator::RFdot(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RB(const auto& k2, const auto& p2)
    {
      return Regulator::dq2RB(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RF(const auto& k2, const auto& p2)
    {
      return Regulator::dq2RF(k2, p2);
    }
  };
} using DiFfRG::Flavour_split_form_kernel;