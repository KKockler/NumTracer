#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class Zq_kernel
  {
    public:
    using Regulator = REG;

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zq_qcd::nenv) > 0 ? (DiFfRG::zq_qcd::nenv) : 1];
      DiFfRG::zq_qcd::fill(fenv, l1, cos1, p);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix1 = ZA.index(l1);
      const auto _interp3 = ZA.at(_ix1);
      const auto _interp4 = ZAAqbq2(0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
      const auto _interp5 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp6 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _ix2 = ZA.index(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp8 = ZA.at(_ix2);
      const auto _ix0 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZAqbq1.at(_ix0);
      const auto _interp10 = zq.at(_ix1);
      const auto _interp11 = ZAqbq4.at(_ix0);
      const auto _interp12 = ZAqbq7.at(_ix0);
      const auto _interp13 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp14 = zq.at(_ix2);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp10 + _interp5 * powr<-1>(l1));
      const auto _den3 = powr<-2>(_interp3 * powr<2>(l1) + _interp1 * nthk0);
      const auto _den4 = powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _den5 = powr<-1>(_interp7 * nthk0 + _interp8 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-1>(_interp14 + _interp13 * sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = _interp2 * nthk0;
      const auto _cse3 = -nthk0;
      const auto _cse4 = _cse3 + nthk2;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + nthk1;
      const auto _cse7 = _cse6 * _interp1;
      const auto _cse8 = _cse2 + _cse7;
      const auto _cse9 = powr<-5>(l1);
      const auto _cse10 = -_interp5 * etaQ;
      const auto _cse11 = _cse10 + _interp6;
      return 0.02083333333333333 * fma(4., _cse1 * _cse11 * _cse9 * _den2 * _den5 * DiFfRG::zq_qcd::tr2(fenv) * powr<2>(_interp11), fma(-4., _cse1 * _cse8 * _den3 * _den4 * _den6 * DiFfRG::zq_qcd::tr6(fenv) * powr<2>(_interp11), fma(-1., _cse1 * _cse11 * _cse9 * _den2 * _den5 * DiFfRG::zq_qcd::tr4(fenv) * powr<2>(_interp12), fma(-1., _cse1 * _cse8 * _den3 * _den4 * _den6 * DiFfRG::zq_qcd::tr8(fenv) * powr<2>(_interp12), fma(2., _cse1 * _cse8 * _den3 * DiFfRG::zq_qcd::tr0(fenv) * _interp4, fma(-2., _cse1 * _cse11 * _cse9 * _den2 * _den5 * DiFfRG::zq_qcd::tr3(fenv) * _interp12 * _interp9, fma(-2., _cse1 * _cse8 * _den3 * _den4 * _den6 * DiFfRG::zq_qcd::tr7(fenv) * _interp12 * _interp9, fma(-4., _cse1 * _cse11 * _cse9 * _den2 * _den5 * DiFfRG::zq_qcd::tr1(fenv) * powr<2>(_interp9), fma(-4., _cse1 * _cse8 * _den3 * _den4 * _den6 * DiFfRG::zq_qcd::tr5(fenv) * powr<2>(_interp9), 0.)))))))));
    }

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const auto _interp1 = zq(p);
      return _interp1 * etaQ;
    }

    // clang-format off
    static device::array<double, 3> ntHoisted(const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      return {{ZA.CPU()(pow(1. + powr<6>(k),0.16666666666666666667)),
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
using DiFfRG::Zq_kernel;