#pragma once

#include "shim.hpp"

namespace DiFfRG {
  template<typename _Regulator>
  class Glu_quark_form_kernel
  {
    public:
    using Regulator = _Regulator;

    static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;const auto _interp1 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp3 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _den1 = powr<-2>(_interp1 + powr<2>(l1));
      const auto _den2 = powr<-1>(_interp2 + powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      return -2.666666666666667 * fma(_den1, _den2 * _interp3 * powr<4>(l1), fma(2., _den1 * _den2 * _interp3 * powr<2>(cos1) * powr<4>(l1), fma(-3., _den1 * _den2 * _interp3 * cos1 * powr<3>(l1) * p, 0.)));
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k)
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
} using DiFfRG::Glu_quark_form_kernel;