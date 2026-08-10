#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZA_kernel
  {
    public:
    using Regulator = REG;

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za_qcd::nenv) > 0 ? (DiFfRG::za_qcd::nenv) : 1];
      DiFfRG::za_qcd::fill(fenv, l1, cos1, p);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix2 = ZA.index(l1);
      const auto _interp3 = ZA.at(_ix2);
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _ix3 = ZA.index(sqrt(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _interp5 = ZA.at(_ix3);
      const auto _ix1 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cos1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA3.at(_ix1);
      const auto _ix0 = ZA4.index(0.7071067811865475 * sqrt(powr<2>(l1) + powr<2>(p)));
      const auto _interp7 = ZA4.at(_ix0);
      const auto _interp8 = ZAcbc.at(_ix1);
      const auto _interp9 = Zc.at(_ix2);
      const auto _interp10 = Zc.at(_ix3);
      const auto _interp11 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp12 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp13 = ZAAqbq2.at(_ix0);
      const auto _interp14 = zq.at(_ix2);
      const auto _interp15 = ZAqbq1.at(_ix1);
      const auto _interp16 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp17 = zq.at(_ix3);
      const auto _interp18 = ZAqbq4.at(_ix1);
      const auto _interp19 = ZAqbq7.at(_ix1);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp14 + _interp11 * powr<-1>(l1));
      const auto _den3 = powr<-2>(_interp3 * powr<2>(l1) + _interp1 * nthk0);
      const auto _den4 = powr<-2>(_interp9 * powr<2>(l1) + _interp1 * nthk3);
      const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _den6 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _den7 = powr<-1>(_interp17 + _interp16 * sqrt(powr<-1>(powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p))));
      const auto _den8 = -powr<-1>(_interp4 * nthk3 + _interp10 * (powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = _interp2 * nthk0;
      const auto _cse3 = -nthk0;
      const auto _cse4 = _cse3 + nthk2;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + nthk1;
      const auto _cse7 = _cse6 * _interp1;
      const auto _cse8 = _cse2 + _cse7;
      const auto _cse9 = powr<-5>(l1);
      const auto _cse10 = -_interp11 * etaQ;
      const auto _cse11 = _cse10 + _interp12;
      return fma(-0.0833333333333333, _cse1 * _cse11 * _cse9 * _den2 * DiFfRG::za_qcd::tr3(fenv) * _interp13, fma(-0.1666666666666667, _cse1 * _cse11 * _cse9 * _den2 * _den5 * _den7 * DiFfRG::za_qcd::tr4(fenv) * powr<2>(_interp15), fma(-0.1666666666666667, _cse1 * _cse11 * _cse9 * _den2 * _den5 * _den7 * DiFfRG::za_qcd::tr5(fenv) * powr<2>(_interp18), fma(-0.0833333333333333, _cse1 * _cse11 * _cse9 * _den2 * _den5 * _den7 * DiFfRG::za_qcd::tr6(fenv) * _interp15 * _interp19, fma(-0.04166666666666666, _cse1 * _cse11 * _cse9 * _den2 * _den5 * _den7 * DiFfRG::za_qcd::tr7(fenv) * powr<2>(_interp19), fma(-0.04166666666666666, _cse1 * _cse8 * _den3 * _den6 * DiFfRG::za_qcd::tr0(fenv) * powr<2>(_interp6), fma(-0.02083333333333333, _cse1 * _cse8 * _den3 * DiFfRG::za_qcd::tr1(fenv) * _interp7, fma(-0.0833333333333333, _cse1 * _den4 * _den8 * DiFfRG::za_qcd::tr2(fenv) * powr<2>(_interp8) * (_interp2 * nthk3 + _interp1 * (nthk4 + 50. * (-nthk3 + nthk5))), 0.))))))));
    }

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    // clang-format on
    {
      return 0.;
    }

    // clang-format off
    static device::array<double, 6> ntHoisted(const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      return {{ZA.CPU()(pow(1. + powr<6>(k),0.16666666666666666667)),
          dtZA.CPU()(pow(1. + powr<6>(k),0.16666666666666666667)),
          ZA.CPU()(1.02 * pow(1. + powr<6>(k),0.16666666666666666667)),
          Zc.CPU()(k),
          dtZc.CPU()(k),
          Zc.CPU()(1.02 * k)}};
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
using DiFfRG::ZA_kernel;