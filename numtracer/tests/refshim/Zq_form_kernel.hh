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
      using namespace DiFfRG;using namespace DiFfRG::compute;const auto _interp1 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp4 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp8 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp9 = Zq(k);
      const auto _interp10 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp11 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cos1 * l1 * p + powr<2>(p)));
      const auto _interp13 = dtZq(k);
      const auto _interp14 = RF(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp15 = RFdot(powr<2>(k), powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp16 = Zq(1.02 * k);
      const auto _interp17 = Mq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = Zq(sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp4 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-1>(_interp2 * _interp4 + _interp6 * powr<2>(l1));
      const auto _den4 = powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _den5 = powr<-1>(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p));
      const auto _den6 = powr<-1>(powr<2>(_interp11) + powr<2>(_interp10) * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)) + _interp8 * _interp9 * (_interp8 * _interp9 + 2. * _interp10 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _den7 = powr<-2>(0.5 * powr<2>(_interp17) + powr<2>(_interp18) * (0.5 * powr<2>(l1) + 1. * cos1 * l1 * p + 0.5 * powr<2>(p)) + _interp14 * _interp9 * (0.5 * _interp14 * _interp9 + 1. * _interp18 * sqrt(powr<2>(l1) + 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = 2. * cos1 * l1 * p;
      const auto _cse2 = _cse1 + powr<2>(l1) + powr<2>(p);
      const auto _cse3 = sqrt(_cse2);
      const auto _cse4 = powr<-1>(p);
      const auto _cse5 = -2. * cos1 * l1 * p;
      const auto _cse6 = _cse5 + powr<2>(l1) + powr<2>(p);
      const auto _cse7 = sqrt(_cse6);
      return fma(1.333333333333333, _cse4 * _cse7 * _den1 * _den2 * _den4 * _den6 * powr<2>(_interp7) * (_cse7 * _interp10 + _interp8 * _interp9) * (_interp2 * (-50. * _interp4 + 50. * _interp5) * powr<6>(k) + _interp3 * _interp4 * (1. + powr<6>(k)) + _interp1 * _interp2 * (1. + 1. * powr<6>(k))) * (-3. * cos1 * l1 + p + 2. * powr<2>(cos1) * p), fma(2., _cse3 * _cse4 * _den3 * _den5 * _den7 * powr<2>(_interp12) * (_interp13 * _interp14 + _interp14 * (50. * _interp16 - 50. * _interp9) + _interp15 * _interp9) * (1. * cos1 * l1 + 0.3333333333333333 * p + 0.6666666666666666 * powr<2>(cos1) * p) * (-0.5 * powr<2>(_interp17) + 1. * _cse3 * _interp14 * _interp18 * _interp9 + 0.5 * powr<2>(_interp14) * powr<2>(_interp9) + powr<2>(_interp18) * (0.5 * powr<2>(l1) + 1. * cos1 * l1 * p + 0.5 * powr<2>(p))), 0.));
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq)
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
} using DiFfRG::Zq_form_kernel;