#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "YM_ZA4_num_kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class YM_ZA4_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::ym_za4_num::nenv) > 0 ? (DiFfRG::ym_za4_num::nenv) : 1];
      DiFfRG::ym_za4_num::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp8 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp10 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp11 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp12 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp14 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp15 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp16 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp17 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp20 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp23 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp24 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp25 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp26 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp27 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp28 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp29 = Zc(k);
      const auto _interp30 = dtZc(k);
      const auto _interp31 = Zc(1.02 * k);
      const auto _interp32 = Zc(l1);
      const auto _interp33 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp34 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp35 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp2 * _interp3 + _interp6 * powr<2>(l1));
      const auto _den3 = powr<-2>(_interp29 * _interp3 + _interp32 * powr<2>(l1));
      const auto _den4 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-1>(_interp2 * _interp21 + _interp22 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den6 = -powr<-1>(_interp29 * _interp7 + _interp33 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den7 = 3. * powr<-1>(3. * _interp11 * _interp2 + _interp12 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den8 = 3. * powr<-1>(3. * _interp18 * _interp2 + _interp19 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den9 = powr<-1>(_interp2 * _interp9 + _interp10 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _den10 = -3. * powr<-1>(3. * _interp11 * _interp29 + _interp35 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den11 = -powr<-1>(_interp29 * _interp9 + _interp34 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _cse1 = _interp1 * _interp2;
      const auto _cse2 = -_interp2;
      const auto _cse3 = _cse2 + _interp5;
      const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
      const auto _cse5 = _cse4 + _interp4;
      const auto _cse6 = _cse5 * _interp3;
      const auto _cse7 = _cse1 + _cse6;
      return fma(0.00510204081632653, _cse7 * _den2 * _den4 * _den7 * _den9 * DiFfRG::ym_za4_num::tr0(fenv) * _interp13 * _interp14 * _interp15 * _interp16, fma(0.00510204081632653, _cse7 * _den2 * _den7 * _den9 * DiFfRG::ym_za4_num::tr1(fenv) * _interp14 * _interp16 * _interp17, fma(0.001275510204081632, _cse7 * _den2 * _den8 * DiFfRG::ym_za4_num::tr3(fenv) * powr<2>(_interp20), fma(0.002551020408163265, _cse7 * _den2 * _den5 * _den9 * DiFfRG::ym_za4_num::tr2(fenv) * _interp14 * _interp23 * _interp24, fma(0.01020408163265306, _den10 * _den11 * _den3 * _den6 * DiFfRG::ym_za4_num::tr4(fenv) * _interp25 * _interp26 * _interp27 * _interp28 * (_interp1 * _interp29 + _interp3 * (_interp30 + 50. * (-_interp29 + _interp31))), 0.)))));
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
using DiFfRG::YM_ZA4_num_kernel;