#pragma once

#include "shim.hpp"

namespace DiFfRG {
  template<typename _Regulator>
  class ZA_form_kernel
  {
    public:
    using Regulator = _Regulator;

    // clang-format off
    static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;
      const auto _interp1 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp4 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
      const auto _interp11 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp12 = dtZc(k);
      const auto _interp13 = Zc(k);
      const auto _interp14 = Zc(1.02 * k);
      const auto _interp15 = Zc(l1);
      const auto _interp16 = Zc(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp18 = dtZq(k);
      const auto _interp19 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp20 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp21 = Zq(k);
      const auto _interp22 = Zq(1.02 * k);
      const auto _interp23 = Mq(l1);
      const auto _interp24 = Zq(l1);
      const auto _interp25 = Mq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp26 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp27 = Zq(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-1>(1. + powr<6>(k));
      const auto _den3 = powr<-2>(_interp13 * _interp2 + _interp15 * powr<2>(l1));
      const auto _den4 = powr<-2>(_interp2 * _interp4 + _interp6 * powr<2>(l1));
      const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _den6 = powr<-1>(_interp13 * _interp7 + _interp16 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den7 = powr<-1>(_interp4 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den8 = powr<-2>(powr<2>(_interp19) * powr<2>(_interp21) + powr<2>(_interp23) + 2. * _interp19 * _interp21 * _interp24 * l1 + powr<2>(_interp24) * powr<2>(l1));
      const auto _den9 = powr<-1>(powr<2>(_interp25) + powr<2>(_interp27) * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)) + _interp21 * _interp26 * (_interp21 * _interp26 + 2. * _interp27 * sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = -2. * cos1 * l1 * p;
      const auto _cse12 = powr<-2>(p);
      const auto _cse2 = _cse1 + powr<-1>(_cse12) + powr<2>(l1);
      const auto _cse3 = sqrt(_cse2);
      const auto _cse4 = -2. * powr<2>(cos1);
      const auto _cse5 = -1. + _cse4;
      const auto _cse6 = _cse5 * l1;
      const auto _cse7 = 3. * cos1 * p;
      const auto _cse8 = _cse6 + _cse7;
      const auto _cse9 = _interp21 * _interp26;
      const auto _cse10 = _cse3 * _interp27;
      const auto _cse11 = _cse10 + _cse9;
      const auto _cse13 = -6. * powr<2>(cos1);
      const auto _cse14 = 6. + _cse13;
      const auto _cse15 = 3. * powr<2>(cos1);
      const auto _cse16 = -3. + _cse15;
      // clang-format off
      return fma(-1., _cse12 * _den4 * _interp10 * (7. - powr<2>(cos1)) * (_interp3 * _interp4 + _interp2 * (_interp1 + 50. * _den1 * (-_interp4 + _interp5) * powr<6>(k))), fma(2., _cse12 * _den3 * _den6 * powr<2>(_interp11) * (_interp12 * _interp2 + (-50. * _interp13 + 50. * _interp14) * _interp2 + _interp13 * _interp3) * (-1. + powr<2>(cos1)) * powr<2>(l1), fma(-4., _cse12 * _den2 * _den4 * _den5 * _den7 * powr<2>(_interp9) * (_interp2 * (-50. * _interp4 + 50. * _interp5) * powr<6>(k) + _interp3 * _interp4 * (1. + powr<6>(k)) + _interp1 * _interp2 * (1. + 1. * powr<6>(k))) * (powr<-2>(_cse12) * _cse16 + powr<-1>(_cse12) * (-8. + 7. * powr<2>(cos1) + powr<4>(cos1)) * powr<2>(l1) + _cse16 * powr<4>(l1) + _cse14 * cos1 * powr<3>(l1) * p + _cse14 * cos1 * l1 * powr<3>(p)), fma(2.666666666666667, _cse12 * _cse3 * _den5 * _den8 * _den9 * powr<2>(_interp17) * (_interp18 * _interp19 + _interp20 * _interp21 + _interp19 * (-50. * _interp21 + 50. * _interp22)) * (-6. * _cse3 * _interp23 * _interp25 * (1. * _interp19 * _interp21 + 1. * _interp24 * l1) + _cse11 * powr<2>(_interp23) * (l1 + 2. * powr<2>(cos1) * l1 - 3. * cos1 * p) + _cse11 * (_cse8 * powr<2>(_interp19) * powr<2>(_interp21) + _cse8 * powr<2>(_interp24) * powr<2>(l1) + _interp19 * _interp21 * _interp24 * l1 * ((-2. - 4. * powr<2>(cos1)) * l1 + 6. * cos1 * p))), 0.))));
      // clang-format on
    }

    // clang-format off
    static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
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
} using DiFfRG::ZA_form_kernel;