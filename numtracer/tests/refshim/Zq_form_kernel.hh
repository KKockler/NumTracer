#pragma once

#include "shim.hpp"

namespace DiFfRG {
  template<typename _Regulator>
  class Zq_form_kernel
  {
    public:
    using Regulator = _Regulator;

    static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;const auto _interp1 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp4 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp5 = dtZq(k);
      const auto _interp6 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp7 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp8 = Zq(k);
      const auto _interp9 = Zq(1.02 * k);
      const auto _interp10 = Mq(l1);
      const auto _interp11 = Zq(l1);
      const auto _interp12 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp13 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp14 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp15 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp16 = ZA(l1);
      const auto _interp17 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp18 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp19 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp13 * _interp2 + _interp16 * powr<2>(l1));
      const auto _den3 = powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _den4 = powr<-1>(_interp1 * _interp2 + _interp3 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-2>(powr<2>(_interp10) + powr<2>(_interp6) * powr<2>(_interp8) + 2. * _interp11 * _interp6 * _interp8 * l1 + powr<2>(_interp11) * powr<2>(l1));
      const auto _den6 = powr<-1>(powr<2>(_interp19) + powr<2>(_interp18) * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)) + _interp17 * _interp8 * (_interp17 * _interp8 + 2. * _interp18 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = -2. * cos1 * l1 * p;
      const auto _cse2 = _cse1 + powr<2>(l1) + powr<2>(p);
      const auto _cse3 = sqrt(_cse2);
      const auto _cse4 = powr<-1>(p);
      return fma(1.333333333333333, _cse3 * _cse4 * _den1 * _den2 * _den3 * _den6 * powr<2>(_interp4) * (_cse3 * _interp18 + _interp17 * _interp8) * (_interp13 * (50. * _interp15 - 50. * _interp2) * powr<6>(k) + _interp14 * _interp2 * (1. + powr<6>(k)) + _interp12 * _interp13 * (1. + 1. * powr<6>(k))) * (-3. * cos1 * l1 + p + 2. * powr<2>(cos1) * p), fma(1.333333333333333, _cse4 * _den3 * _den4 * _den5 * powr<2>(_interp4) * (_interp5 * _interp6 + _interp7 * _interp8 + _interp6 * (-50. * _interp8 + 50. * _interp9)) * (-powr<2>(_interp10) + powr<2>(_interp6) * powr<2>(_interp8) + 2. * _interp11 * _interp6 * _interp8 * l1 + powr<2>(_interp11) * powr<2>(l1)) * (3. * cos1 * powr<2>(l1) - 2. * l1 * p - 4. * powr<2>(cos1) * l1 * p + 3. * cos1 * powr<2>(p)), 0.));
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;
      return 0.;
    }
    private:
    static KOKKOS_FORCEINLINE_FUNCTION auto RB(const auto& k2, const auto& p2)
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
} using DiFfRG::Zq_form_kernel;