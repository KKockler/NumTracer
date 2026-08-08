#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "YM_ZA4_num_kernels.hh"

namespace DiFfRG
{
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
      const auto _interp7 = RB(powr<2>(k), fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp8 = ZA(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp9 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp10 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp11 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp12 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp14 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp15 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp16 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp17 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp18 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp19 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp20 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp21 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp22 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp23 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp24 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp25 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp26 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp28 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp29 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp30 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp31 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp32 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp33 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp34 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp35 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp36 = ZA3(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp38 = RB(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp39 = ZA(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp40 = ZA3(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp42 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp43 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp44 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp46 = RB(powr<2>(k), fma(2., cosl1p1 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp47 = ZA(sqrt(fma(2., cosl1p1 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp48 = ZA3(0.816496580927726 * sqrt(fma(cosl1p1, l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp62 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp63 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp66 = Zc(k);
      const auto _interp67 = dtZc(k);
      const auto _interp68 = Zc(1.02 * k);
      const auto _interp69 = Zc(l1);
      const auto _interp70 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp71 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp73 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp75 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp76 = ZAcbc(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp78 = Zc(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp79 = ZAcbc(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp82 = Zc(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp85 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp87 = Zc(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp88 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp89 = ZAcbc(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp90 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp2, _interp3, fma(_interp6, powr<2>(l1), 0.)));
      const auto _den3 = powr<-2>(fma(-1., _interp3 * _interp66, fma(-1., _interp69 * powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp2, _interp21, fma(_interp22, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp2, _interp25, fma(_interp26, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den6 = powr<-1>(fma(_interp2, _interp34, fma(_interp35, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den7 = powr<-1>(fma(_interp2, _interp7, fma(_interp8, powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(_interp2, _interp38, fma(_interp39, powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den9 = powr<-1>(fma(-1., _interp21 * _interp66, fma(-1., _interp70 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)), 0.)));
      const auto _den10 = powr<-1>(fma(-1., _interp25 * _interp66, fma(-1., _interp75 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)), 0.)));
      const auto _den11 = powr<-1>(fma(-1., _interp34 * _interp66, fma(-1., _interp82 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)), 0.)));
      const auto _den12 = powr<-1>(fma(-1., _interp66 * _interp7, fma(-1., _interp87 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)), 0.)));
      const auto _den13 = powr<-1>(fma(-1., _interp38 * _interp66, fma(-1., _interp78 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)), 0.)));
      const auto _den14 = powr<-1>(fma(_interp2, _interp29, fma(_interp30, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den15 = powr<-1>(fma(_interp17, _interp2, fma(_interp18, powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den16 = powr<-1>(fma(_interp11, _interp2, fma(_interp12, powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den17 = powr<-1>(fma(_interp2, _interp46, fma(_interp47, powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den18 = powr<-1>(fma(_interp2, _interp42, fma(_interp43, powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den19 = powr<-1>(fma(_interp2, _interp9, fma(_interp10, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
      const auto _den22 = powr<-1>(fma(-1., _interp11 * _interp66, fma(-1., _interp90 * powr<2>(l1) - 2. * _interp90 * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * _interp90 * powr<2>(p), 0.)));
      const auto _den25 = powr<-1>(fma(-1., _interp66 * _interp9, fma(-1., _interp71 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)), 0.)));
      const auto _interp64 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp65 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp72 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp74 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp77 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp80 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp81 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp83 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp84 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp86 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _den20 = -3. * powr<-1>(3. * _interp29 * _interp66 + _interp72 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den21 = -3. * powr<-1>(3. * _interp17 * _interp66 + _interp83 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _interp91 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp92 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp93 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp94 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp95 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp96 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp97 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp98 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _den23 = -3. * powr<-1>(3. * _interp46 * _interp66 + _interp94 * (3. * powr<2>(l1) + 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den24 = -3. * powr<-1>(3. * _interp42 * _interp66 + _interp97 * (3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _interp37 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp41 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp45 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp49 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp50 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp51 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp52 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp53 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp54 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p1 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp55 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp56 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp57 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp58 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp59 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (2. * cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp60 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp61 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den11 + _den12 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den2 + _den20 + _den21 + _den22 + _den23 + _den24 + _den25 + _den3 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + DiFfRG::ym_za4_num::tr0(fenv) + DiFfRG::ym_za4_num::tr10(fenv) + DiFfRG::ym_za4_num::tr11(fenv) + DiFfRG::ym_za4_num::tr12(fenv) + DiFfRG::ym_za4_num::tr13(fenv) + DiFfRG::ym_za4_num::tr14(fenv) + DiFfRG::ym_za4_num::tr15(fenv) + DiFfRG::ym_za4_num::tr16(fenv) + DiFfRG::ym_za4_num::tr17(fenv) + DiFfRG::ym_za4_num::tr18(fenv) + DiFfRG::ym_za4_num::tr19(fenv) + DiFfRG::ym_za4_num::tr1(fenv) + DiFfRG::ym_za4_num::tr20(fenv) + DiFfRG::ym_za4_num::tr21(fenv) + DiFfRG::ym_za4_num::tr22(fenv) + DiFfRG::ym_za4_num::tr23(fenv) + DiFfRG::ym_za4_num::tr24(fenv) + DiFfRG::ym_za4_num::tr25(fenv) + DiFfRG::ym_za4_num::tr26(fenv) + DiFfRG::ym_za4_num::tr27(fenv) + DiFfRG::ym_za4_num::tr28(fenv) + DiFfRG::ym_za4_num::tr29(fenv) + DiFfRG::ym_za4_num::tr2(fenv) + DiFfRG::ym_za4_num::tr30(fenv) + DiFfRG::ym_za4_num::tr31(fenv) + DiFfRG::ym_za4_num::tr32(fenv) + DiFfRG::ym_za4_num::tr33(fenv) + DiFfRG::ym_za4_num::tr34(fenv) + DiFfRG::ym_za4_num::tr35(fenv) + DiFfRG::ym_za4_num::tr36(fenv) + DiFfRG::ym_za4_num::tr37(fenv) + DiFfRG::ym_za4_num::tr38(fenv) + DiFfRG::ym_za4_num::tr39(fenv) + DiFfRG::ym_za4_num::tr3(fenv) + DiFfRG::ym_za4_num::tr40(fenv) + DiFfRG::ym_za4_num::tr41(fenv) + DiFfRG::ym_za4_num::tr42(fenv) + DiFfRG::ym_za4_num::tr43(fenv) + DiFfRG::ym_za4_num::tr44(fenv) + DiFfRG::ym_za4_num::tr45(fenv) + DiFfRG::ym_za4_num::tr46(fenv) + DiFfRG::ym_za4_num::tr47(fenv) + DiFfRG::ym_za4_num::tr48(fenv) + DiFfRG::ym_za4_num::tr49(fenv) + DiFfRG::ym_za4_num::tr4(fenv) + DiFfRG::ym_za4_num::tr50(fenv) + DiFfRG::ym_za4_num::tr51(fenv) + DiFfRG::ym_za4_num::tr52(fenv) + DiFfRG::ym_za4_num::tr53(fenv) + DiFfRG::ym_za4_num::tr54(fenv) + DiFfRG::ym_za4_num::tr55(fenv) + DiFfRG::ym_za4_num::tr56(fenv) + DiFfRG::ym_za4_num::tr5(fenv) + DiFfRG::ym_za4_num::tr6(fenv) + DiFfRG::ym_za4_num::tr7(fenv) + DiFfRG::ym_za4_num::tr8(fenv) + DiFfRG::ym_za4_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + k);
      // clang-format on
      _T _acc{};
      { // subkernel 1
      // clang-format off
      _acc += fma(-0.0004251700680272108, _den20 * _den25 * _den3 * _den9 * DiFfRG::ym_za4_num::tr45(fenv) * _interp62 * _interp63 * _interp64 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))), fma(0.0004251700680272108, _den20 * _den25 * _den3 * _den9 * DiFfRG::ym_za4_num::tr46(fenv) * _interp62 * _interp63 * _interp64 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))), fma(-0.0004251700680272108, _den10 * _den20 * _den25 * _den3 * DiFfRG::ym_za4_num::tr47(fenv) * _interp63 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74, fma(0.0004251700680272108, _den10 * _den20 * _den25 * _den3 * DiFfRG::ym_za4_num::tr48(fenv) * _interp63 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74, fma(0.0004251700680272108, _den13 * _den20 * _den3 * _den9 * DiFfRG::ym_za4_num::tr33(fenv) * _interp62 * _interp64 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp77, fma(-0.0004251700680272108, _den13 * _den20 * _den3 * _den9 * DiFfRG::ym_za4_num::tr34(fenv) * _interp62 * _interp64 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp77, fma(0.0004251700680272108, _den10 * _den13 * _den20 * _den3 * DiFfRG::ym_za4_num::tr35(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74 * _interp76 * _interp77, fma(-0.0004251700680272108, _den10 * _den13 * _den20 * _den3 * DiFfRG::ym_za4_num::tr36(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74 * _interp76 * _interp77, fma(-0.0004251700680272108, _den11 * _den21 * _den3 * _den9 * DiFfRG::ym_za4_num::tr37(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp80 * _interp81, fma(0.0004251700680272108, _den11 * _den21 * _den3 * _den9 * DiFfRG::ym_za4_num::tr38(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp80 * _interp81, fma(-0.0004251700680272108, _den21 * _den25 * _den3 * _den9 * DiFfRG::ym_za4_num::tr49(fenv) * _interp62 * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp80 * _interp84, fma(0.0004251700680272108, _den21 * _den25 * _den3 * _den9 * DiFfRG::ym_za4_num::tr50(fenv) * _interp62 * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp80 * _interp84, fma(-0.0004251700680272108, _den12 * _den21 * _den25 * _den3 * DiFfRG::ym_za4_num::tr51(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp84 * _interp85 * _interp86, fma(0.0004251700680272108, _den12 * _den21 * _den25 * _den3 * DiFfRG::ym_za4_num::tr52(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp84 * _interp85 * _interp86, fma(-0.0004251700680272108, _den10 * _den22 * _den25 * _den3 * DiFfRG::ym_za4_num::tr53(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp88 * _interp89, 0.)))))))))))))));
      // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(0.0004251700680272108, _den10 * _den22 * _den25 * _den3 * DiFfRG::ym_za4_num::tr54(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp88 * _interp89, fma(-0.0004251700680272108, _den12 * _den22 * _den25 * _den3 * DiFfRG::ym_za4_num::tr55(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp85 * _interp88 * _interp91, fma(0.0004251700680272108, _den12 * _den22 * _den25 * _den3 * DiFfRG::ym_za4_num::tr56(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp85 * _interp88 * _interp91, fma(0.0004251700680272108, _den10 * _den13 * _den23 * _den3 * DiFfRG::ym_za4_num::tr39(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp76 * _interp92 * _interp93, fma(-0.0004251700680272108, _den10 * _den13 * _den23 * _den3 * DiFfRG::ym_za4_num::tr40(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp76 * _interp92 * _interp93, fma(0.0004251700680272108, _den11 * _den24 * _den3 * _den9 * DiFfRG::ym_za4_num::tr41(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp95 * _interp96, fma(-0.0004251700680272108, _den11 * _den24 * _den3 * _den9 * DiFfRG::ym_za4_num::tr42(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp95 * _interp96, fma(0.0004251700680272108, _den13 * _den24 * _den3 * _den9 * DiFfRG::ym_za4_num::tr43(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp95 * _interp98, fma(-0.0004251700680272108, _den13 * _den24 * _den3 * _den9 * DiFfRG::ym_za4_num::tr44(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp95 * _interp98, fma(0.0004251700680272108, _den16 * _den19 * _den2 * _den7 * DiFfRG::ym_za4_num::tr10(fenv) * _interp13 * _interp14 * _interp15 * _interp16 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den19 * _den2 * _den7 * DiFfRG::ym_za4_num::tr8(fenv) * _interp13 * _interp14 * _interp19 * _interp20 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den19 * _den2 * _den4 * DiFfRG::ym_za4_num::tr9(fenv) * _interp14 * _interp20 * _interp23 * _interp24 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den19 * _den2 * _den5 * DiFfRG::ym_za4_num::tr11(fenv) * _interp14 * _interp15 * _interp27 * _interp28 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den19 * _den2 * _den5 * DiFfRG::ym_za4_num::tr6(fenv) * _interp14 * _interp27 * _interp31 * _interp32 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den19 * _den2 * _den4 * DiFfRG::ym_za4_num::tr7(fenv) * _interp14 * _interp23 * _interp32 * _interp33 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), 0.)))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.0004251700680272108, _den15 * _den2 * _den4 * _den6 * DiFfRG::ym_za4_num::tr2(fenv) * _interp23 * _interp24 * _interp36 * _interp37 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den5 * _den8 * DiFfRG::ym_za4_num::tr0(fenv) * _interp27 * _interp31 * _interp40 * _interp41 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den4 * _den8 * DiFfRG::ym_za4_num::tr1(fenv) * _interp23 * _interp33 * _interp40 * _interp41 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den4 * _den6 * DiFfRG::ym_za4_num::tr4(fenv) * _interp23 * _interp36 * _interp44 * _interp45 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den5 * _den8 * DiFfRG::ym_za4_num::tr3(fenv) * _interp27 * _interp40 * _interp48 * _interp49 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den4 * _den8 * DiFfRG::ym_za4_num::tr5(fenv) * _interp23 * _interp40 * _interp44 * _interp50 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den5 * DiFfRG::ym_za4_num::tr24(fenv) * _interp27 * _interp31 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den19 * _den2 * DiFfRG::ym_za4_num::tr13(fenv) * _interp14 * _interp32 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den4 * DiFfRG::ym_za4_num::tr25(fenv) * _interp23 * _interp33 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den8 * DiFfRG::ym_za4_num::tr12(fenv) * _interp40 * _interp41 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * DiFfRG::ym_za4_num::tr26(fenv) * powr<2>(_interp51) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den7 * DiFfRG::ym_za4_num::tr27(fenv) * _interp13 * _interp19 * _interp52 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den19 * _den2 * DiFfRG::ym_za4_num::tr14(fenv) * _interp14 * _interp20 * _interp52 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den4 * DiFfRG::ym_za4_num::tr28(fenv) * _interp23 * _interp24 * _interp52 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * DiFfRG::ym_za4_num::tr29(fenv) * powr<2>(_interp52) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), 0.)))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(0.0004251700680272108, _den16 * _den19 * _den2 * DiFfRG::ym_za4_num::tr16(fenv) * _interp14 * _interp15 * _interp53 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den7 * DiFfRG::ym_za4_num::tr30(fenv) * _interp13 * _interp16 * _interp53 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den5 * DiFfRG::ym_za4_num::tr31(fenv) * _interp27 * _interp28 * _interp53 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * DiFfRG::ym_za4_num::tr32(fenv) * powr<2>(_interp53) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den5 * DiFfRG::ym_za4_num::tr15(fenv) * _interp27 * _interp48 * _interp54 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den4 * DiFfRG::ym_za4_num::tr17(fenv) * _interp23 * _interp44 * _interp55 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den6 * DiFfRG::ym_za4_num::tr18(fenv) * _interp23 * _interp36 * _interp56 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den7 * DiFfRG::ym_za4_num::tr19(fenv) * _interp13 * _interp14 * _interp57 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den5 * DiFfRG::ym_za4_num::tr20(fenv) * _interp14 * _interp27 * _interp58 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * DiFfRG::ym_za4_num::tr21(fenv) * _interp14 * _interp23 * _interp59 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den8 * DiFfRG::ym_za4_num::tr22(fenv) * _interp23 * _interp40 * _interp60 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den5 * _den8 * DiFfRG::ym_za4_num::tr23(fenv) * _interp27 * _interp40 * _interp61 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), 0.))))))))))));
      // clang-format on
      }
      return _acc;
    }

    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      return 0.;
    }
  };
}
using DiFfRG::YM_ZA4_num_kernel;