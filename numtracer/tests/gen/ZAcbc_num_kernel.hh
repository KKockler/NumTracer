#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "ZAcbc_num_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZAcbc_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zacbc_num::nenv) > 0 ? (DiFfRG::zacbc_num::nenv) : 1];
      DiFfRG::zacbc_num::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp5 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp7 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp8 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp9 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp10 = Zc(k);
      const auto _interp11 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp12 = dtZc(k);
      const auto _interp13 = Zc(1.02 * k);
      const auto _interp14 = Zc(l1);
      const auto _interp15 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp16 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp17 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp18 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp19 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp20 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp21 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp22 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp23 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp24 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp25 = ZA(l1);
      const auto _interp26 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp27 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp28 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp29 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp30 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp31 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp11 * _interp2 + _interp25 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp10 * _interp11 + _interp14 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp17 * _interp2 + _interp26 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-1>(_interp1 * _interp2 + _interp3 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-1>(_interp19 * _interp2 + _interp20 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den7 = -powr<-1>(_interp10 * _interp17 + _interp18 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den8 = -powr<-1>(_interp1 * _interp10 + _interp30 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den9 = -powr<-1>(_interp10 * _interp19 + _interp31 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den10 = powr<-1>(_interp2 * _interp4 + _interp5 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den11 = -powr<-1>(_interp10 * _interp4 + _interp28 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = _interp2 * _interp9;
      const auto _cse3 = -_interp2;
      const auto _cse4 = _cse3 + _interp24;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + _interp23;
      const auto _cse7 = _cse6 * _interp11;
      const auto _cse8 = _cse2 + _cse7;
      const auto _cse9 = _interp10 * _interp9;
      const auto _cse10 = -_interp10;
      const auto _cse11 = _cse10 + _interp13;
      const auto _cse12 = 50. * _cse11;
      const auto _cse13 = _cse12 + _interp12;
      const auto _cse14 = _cse13 * _interp11;
      const auto _cse15 = _cse14 + _cse9;
      return 0.05555555555555555 * fma(_cse1, _cse15 * _den3 * _den6 * _den7 * DiFfRG::zacbc_num::tr1(fenv) * _interp15 * _interp21 * _interp22, fma(_cse1, _cse8 * _den2 * _den4 * _den9 * DiFfRG::zacbc_num::tr3(fenv) * _interp21 * _interp22 * _interp27, fma(-1., _cse1 * _cse15 * _den10 * _den3 * _den7 * DiFfRG::zacbc_num::tr2(fenv) * _interp15 * _interp16 * _interp8, fma(_cse1, _cse8 * _den11 * _den2 * _den4 * DiFfRG::zacbc_num::tr4(fenv) * _interp16 * _interp27 * _interp8, fma(_cse1, _cse8 * _den11 * _den2 * _den8 * DiFfRG::zacbc_num::tr5(fenv) * _interp29 * _interp7 * _interp8, fma(-1., _cse1 * _cse15 * _den10 * _den3 * _den5 * DiFfRG::zacbc_num::tr0(fenv) * _interp6 * _interp7 * _interp8, 0.))))));
    }

    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      return 0.;
    }
    private:
    static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using DiFfRG::ZAcbc_num_kernel;