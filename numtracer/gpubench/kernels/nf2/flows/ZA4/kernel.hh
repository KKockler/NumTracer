#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZA4_kernel
  {
    public:
    using Regulator = REG;

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za4_qcd::nenv) > 0 ? (DiFfRG::za4_qcd::nenv) : 1];
      DiFfRG::za4_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp30 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp31 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _ix5 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp34 = ZAqbq1.at(_ix5);
      const auto _ix12 = zq.index(l1);
      const auto _interp35 = zq.at(_ix12);
      const auto _interp38 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp39 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp42 = ZAqbq4.at(_ix5);
      const auto _interp44 = ZAqbq7.at(_ix5);
      const auto _ix3 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp45 = ZAAqbq2.at(_ix3);
      const auto _ix8 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp46 = ZAqbq1.at(_ix8);
      const auto _interp47 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp48 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp49 = ZAAqbq1.at(_ix3);
      const auto _interp50 = ZAqbq4.at(_ix8);
      const auto _interp51 = ZAqbq7.at(_ix8);
      const auto _ix6 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp52 = ZAqbq1.at(_ix6);
      const auto _ix7 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp53 = ZAqbq1.at(_ix7);
      const auto _interp54 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp55 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp56 = ZAqbq4.at(_ix6);
      const auto _interp57 = ZAqbq4.at(_ix7);
      const auto _interp58 = ZAqbq7.at(_ix6);
      const auto _interp59 = ZAqbq7.at(_ix7);
      const auto _den2 = powr<-2>(fma(_interp30, powr<-1>(l1), _interp35));
      const auto _den5 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den9 = powr<-1>(fma(_interp54, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp55));
      const auto _den12 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _den14 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _den18 = powr<-1>(fma(_interp47, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))), _interp48));
      const auto _den21 = powr<-1>(fma(_interp38, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))), _interp39));
      const auto _ix4 = ZAAqbq2.index(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp32 = ZAAqbq2.at(_ix4);
      const auto _ix11 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp33 = ZAqbq1.at(_ix11);
      const auto _interp36 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix17 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp37 = zq.at(_ix17);
      const auto _interp40 = ZAAqbq1.at(_ix4);
      const auto _interp41 = ZAqbq4.at(_ix11);
      const auto _interp43 = ZAqbq7.at(_ix11);
      const auto _den11 = powr<-1>(_interp37 + _interp36 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA.at(_ix12);
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix15 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp5 = ZA.at(_ix15);
      const auto _interp6 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _ix13 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp7 = ZA.at(_ix13);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix14 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp9 = ZA.at(_ix14);
      const auto _ix10 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3.at(_ix10);
      const auto _ix9 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp11 = ZA3.at(_ix9);
      const auto _ix1 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp12 = ZA3.at(_ix1);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp13 = ZA3.at(_ix0);
      const auto _interp14 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix16 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp16 = ZA.at(_ix16);
      const auto _ix2 = ZA4.index(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp17 = ZA4.at(_ix2);
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp19 = ZA.at(_ix17);
      const auto _interp20 = ZA3.at(_ix11);
      const auto _interp21 = ZA4.at(_ix4);
      const auto _interp22 = ZAcbc.at(_ix10);
      const auto _interp23 = ZAcbc.at(_ix9);
      const auto _interp24 = ZAcbc.at(_ix1);
      const auto _interp25 = ZAcbc.at(_ix0);
      const auto _interp26 = Zc.at(_ix12);
      const auto _interp27 = Zc.at(_ix15);
      const auto _interp28 = Zc.at(_ix13);
      const auto _interp29 = Zc.at(_ix14);
      const auto _interp60 = ZAAqbq1.at(_ix2);
      const auto _interp61 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp62 = zq.at(_ix16);
      const auto _interp63 = ZAAqbq2.at(_ix2);
      const auto _den13 = powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den15 = 3. * powr<-1>(3. * _interp8 * nthk0 + _interp9 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den16 = 3. * powr<-1>(3. * _interp15 * nthk0 + _interp16 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den19 = -3. * powr<-1>(3. * _interp8 * nthk3 + _interp29 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den20 = powr<-1>(_interp62 + _interp61 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      // clang-format off
      using _T = decltype(_den11 + _den12 + _den13 + _den14 + _den15 + _den16 + _den18 + _den19 + _den2 + _den20 + _den21 + _den5 + _den9 + DiFfRG::za4_qcd::tr0(fenv) + DiFfRG::za4_qcd::tr10(fenv) + DiFfRG::za4_qcd::tr11(fenv) + DiFfRG::za4_qcd::tr12(fenv) + DiFfRG::za4_qcd::tr13(fenv) + DiFfRG::za4_qcd::tr14(fenv) + DiFfRG::za4_qcd::tr15(fenv) + DiFfRG::za4_qcd::tr16(fenv) + DiFfRG::za4_qcd::tr17(fenv) + DiFfRG::za4_qcd::tr18(fenv) + DiFfRG::za4_qcd::tr19(fenv) + DiFfRG::za4_qcd::tr1(fenv) + DiFfRG::za4_qcd::tr20(fenv) + DiFfRG::za4_qcd::tr21(fenv) + DiFfRG::za4_qcd::tr22(fenv) + DiFfRG::za4_qcd::tr23(fenv) + DiFfRG::za4_qcd::tr24(fenv) + DiFfRG::za4_qcd::tr25(fenv) + DiFfRG::za4_qcd::tr26(fenv) + DiFfRG::za4_qcd::tr27(fenv) + DiFfRG::za4_qcd::tr28(fenv) + DiFfRG::za4_qcd::tr29(fenv) + DiFfRG::za4_qcd::tr2(fenv) + DiFfRG::za4_qcd::tr30(fenv) + DiFfRG::za4_qcd::tr31(fenv) + DiFfRG::za4_qcd::tr32(fenv) + DiFfRG::za4_qcd::tr33(fenv) + DiFfRG::za4_qcd::tr34(fenv) + DiFfRG::za4_qcd::tr35(fenv) + DiFfRG::za4_qcd::tr36(fenv) + DiFfRG::za4_qcd::tr37(fenv) + DiFfRG::za4_qcd::tr38(fenv) + DiFfRG::za4_qcd::tr39(fenv) + DiFfRG::za4_qcd::tr3(fenv) + DiFfRG::za4_qcd::tr40(fenv) + DiFfRG::za4_qcd::tr41(fenv) + DiFfRG::za4_qcd::tr42(fenv) + DiFfRG::za4_qcd::tr43(fenv) + DiFfRG::za4_qcd::tr44(fenv) + DiFfRG::za4_qcd::tr45(fenv) + DiFfRG::za4_qcd::tr46(fenv) + DiFfRG::za4_qcd::tr47(fenv) + DiFfRG::za4_qcd::tr48(fenv) + DiFfRG::za4_qcd::tr49(fenv) + DiFfRG::za4_qcd::tr4(fenv) + DiFfRG::za4_qcd::tr50(fenv) + DiFfRG::za4_qcd::tr51(fenv) + DiFfRG::za4_qcd::tr52(fenv) + DiFfRG::za4_qcd::tr53(fenv) + DiFfRG::za4_qcd::tr54(fenv) + DiFfRG::za4_qcd::tr55(fenv) + DiFfRG::za4_qcd::tr56(fenv) + DiFfRG::za4_qcd::tr57(fenv) + DiFfRG::za4_qcd::tr58(fenv) + DiFfRG::za4_qcd::tr59(fenv) + DiFfRG::za4_qcd::tr5(fenv) + DiFfRG::za4_qcd::tr60(fenv) + DiFfRG::za4_qcd::tr61(fenv) + DiFfRG::za4_qcd::tr62(fenv) + DiFfRG::za4_qcd::tr63(fenv) + DiFfRG::za4_qcd::tr64(fenv) + DiFfRG::za4_qcd::tr65(fenv) + DiFfRG::za4_qcd::tr66(fenv) + DiFfRG::za4_qcd::tr67(fenv) + DiFfRG::za4_qcd::tr68(fenv) + DiFfRG::za4_qcd::tr69(fenv) + DiFfRG::za4_qcd::tr6(fenv) + DiFfRG::za4_qcd::tr70(fenv) + DiFfRG::za4_qcd::tr71(fenv) + DiFfRG::za4_qcd::tr72(fenv) + DiFfRG::za4_qcd::tr73(fenv) + DiFfRG::za4_qcd::tr74(fenv) + DiFfRG::za4_qcd::tr7(fenv) + DiFfRG::za4_qcd::tr8(fenv) + DiFfRG::za4_qcd::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + nthk3 + nthk4 + nthk5 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den6 = powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _cse1_k1 = powr<-5>(l1);
        const auto _cse2_k1 = -_interp30 * etaQ;
        const auto _cse3_k1 = _cse2_k1 + _interp31;
        // clang-format off
        _acc += 0.002551020408163265 * fma(4., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr7(fenv) * _interp32 * _interp33 * _interp34, fma(-4., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr8(fenv) * _interp34 * _interp40 * _interp41, fma(4., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr9(fenv) * _interp33 * _interp40 * _interp42, fma(4., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr10(fenv) * _interp32 * _interp41 * _interp42, fma(2., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr11(fenv) * _interp32 * _interp34 * _interp43, fma(2., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr12(fenv) * _interp40 * _interp42 * _interp43, fma(2., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr13(fenv) * _interp32 * _interp33 * _interp44, fma(-2., _cse1_k1 * _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr14(fenv) * _interp40 * _interp41 * _interp44, fma(_cse1_k1, _cse3_k1 * _den11 * _den14 * _den2 * _den21 * _den6 * DiFfRG::za4_qcd::tr15(fenv) * _interp32 * _interp43 * _interp44, fma(4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr16(fenv) * _interp34 * _interp45 * _interp46, fma(-4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr17(fenv) * _interp34 * _interp45 * _interp46, fma(2., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr28(fenv) * _interp44 * _interp45 * _interp46, fma(-2., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr29(fenv) * _interp44 * _interp45 * _interp46, fma(4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr20(fenv) * _interp42 * _interp46 * _interp49, fma(-4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr21(fenv) * _interp42 * _interp46 * _interp49, fma(4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr22(fenv) * _interp42 * _interp45 * _interp50, fma(-4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr23(fenv) * _interp42 * _interp45 * _interp50, fma(-4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr18(fenv) * _interp34 * _interp49 * _interp50, fma(4., _cse1_k1 * _cse3_k1 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr19(fenv) * _interp34 * _interp49 * _interp50, 0.)))))))))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _cse1_k2 = powr<-5>(l1);
        const auto _cse2_k2 = -_interp30 * etaQ;
        const auto _cse3_k2 = _cse2_k2 + _interp31;
        // clang-format off
        _acc += 0.002551020408163265 * fma(-2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr30(fenv) * _interp44 * _interp49 * _interp50, fma(2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr31(fenv) * _interp44 * _interp49 * _interp50, fma(2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr24(fenv) * _interp34 * _interp45 * _interp51, fma(-2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr25(fenv) * _interp34 * _interp45 * _interp51, fma(_cse1_k2, _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr32(fenv) * _interp44 * _interp45 * _interp51, fma(-1., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr33(fenv) * _interp44 * _interp45 * _interp51, fma(2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr26(fenv) * _interp42 * _interp49 * _interp51, fma(-2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * DiFfRG::za4_qcd::tr27(fenv) * _interp42 * _interp49 * _interp51, fma(-8., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr34(fenv) * _interp34 * _interp46 * _interp52 * _interp53, fma(-4., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr61(fenv) * _interp44 * _interp46 * _interp52 * _interp53, fma(-8., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr40(fenv) * _interp42 * _interp50 * _interp52 * _interp53, fma(-4., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr52(fenv) * _interp34 * _interp51 * _interp52 * _interp53, fma(-2., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr70(fenv) * _interp44 * _interp51 * _interp52 * _interp53, fma(-8., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr38(fenv) * _interp42 * _interp46 * _interp53 * _interp56, fma(8., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr36(fenv) * _interp34 * _interp50 * _interp53 * _interp56, fma(4., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr63(fenv) * _interp44 * _interp50 * _interp53 * _interp56, fma(-4., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr54(fenv) * _interp42 * _interp51 * _interp53 * _interp56, fma(-8., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr39(fenv) * _interp42 * _interp46 * _interp52 * _interp57, fma(8., _cse1_k2 * _cse3_k2 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr37(fenv) * _interp34 * _interp50 * _interp52 * _interp57, 0.)))))))))))))))))));
        // clang-format on
      }
      { // subkernel 3
        const auto _cse1_k3 = powr<-5>(l1);
        const auto _cse2_k3 = -_interp30 * etaQ;
        const auto _cse3_k3 = _cse2_k3 + _interp31;
        // clang-format off
        _acc += 0.002551020408163265 * fma(4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr64(fenv) * _interp44 * _interp50 * _interp52 * _interp57, fma(-4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr55(fenv) * _interp42 * _interp51 * _interp52 * _interp57, fma(8., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr35(fenv) * _interp34 * _interp46 * _interp56 * _interp57, fma(4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr62(fenv) * _interp44 * _interp46 * _interp56 * _interp57, fma(8., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr41(fenv) * _interp42 * _interp50 * _interp56 * _interp57, fma(4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr53(fenv) * _interp34 * _interp51 * _interp56 * _interp57, fma(2., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr71(fenv) * _interp44 * _interp51 * _interp56 * _interp57, fma(-4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr42(fenv) * _interp34 * _interp46 * _interp53 * _interp58, fma(-2., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr65(fenv) * _interp44 * _interp46 * _interp53 * _interp58, fma(-4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr45(fenv) * _interp42 * _interp50 * _interp53 * _interp58, fma(-2., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr56(fenv) * _interp34 * _interp51 * _interp53 * _interp58, fma(-1., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr72(fenv) * _interp44 * _interp51 * _interp53 * _interp58, fma(-4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr44(fenv) * _interp42 * _interp46 * _interp57 * _interp58, fma(4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr43(fenv) * _interp34 * _interp50 * _interp57 * _interp58, fma(2., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr66(fenv) * _interp44 * _interp50 * _interp57 * _interp58, fma(-2., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr57(fenv) * _interp42 * _interp51 * _interp57 * _interp58, fma(-4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr46(fenv) * _interp34 * _interp46 * _interp52 * _interp59, fma(-2., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr67(fenv) * _interp44 * _interp46 * _interp52 * _interp59, fma(-4., _cse1_k3 * _cse3_k3 * _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr49(fenv) * _interp42 * _interp50 * _interp52 * _interp59, 0.)))))))))))))))))));
        // clang-format on
      }
      { // subkernel 4
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den3 = powr<-2>(_interp3 * powr<2>(l1) + _interp1 * nthk0);
        const auto _den4 = powr<-2>(_interp26 * powr<2>(l1) + _interp1 * nthk3);
        const auto _den7 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den8 = powr<-1>(_interp18 * nthk0 + _interp19 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den10 = -powr<-1>(_interp4 * nthk3 + _interp27 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den17 = powr<-1>(_interp6 * nthk0 + _interp7 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        const auto _den22 = -powr<-1>(_interp6 * nthk3 + _interp28 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(-0.00510204081632653, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr58(fenv) * _interp34 * _interp51 * _interp52 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.002551020408163265, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr73(fenv) * _interp44 * _interp51 * _interp52 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.01020408163265306, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr48(fenv) * _interp42 * _interp46 * _interp56 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(0.01020408163265306, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr47(fenv) * _interp34 * _interp50 * _interp56 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(0.00510204081632653, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr68(fenv) * _interp44 * _interp50 * _interp56 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.00510204081632653, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr59(fenv) * _interp42 * _interp51 * _interp56 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.00510204081632653, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr50(fenv) * _interp34 * _interp46 * _interp58 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.002551020408163265, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr69(fenv) * _interp44 * _interp46 * _interp58 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.00510204081632653, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr51(fenv) * _interp42 * _interp50 * _interp58 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.002551020408163265, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr60(fenv) * _interp34 * _interp51 * _interp58 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.001275510204081632, _den12 * _den14 * _den18 * _den2 * _den21 * _den5 * _den9 * DiFfRG::za4_qcd::tr74(fenv) * _interp44 * _interp51 * _interp58 * _interp59 * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(0.00510204081632653, _den13 * _den2 * _den20 * DiFfRG::za4_qcd::tr5(fenv) * powr<2>(_interp60) * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(-0.00510204081632653, _den13 * _den2 * _den20 * DiFfRG::za4_qcd::tr6(fenv) * powr<2>(_interp63) * (_interp31 - _interp30 * etaQ) * powr<-5>(l1), fma(0.00510204081632653, _den15 * _den17 * _den3 * _den7 * DiFfRG::za4_qcd::tr0(fenv) * _interp10 * _interp11 * _interp12 * _interp13 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.00510204081632653, _den15 * _den17 * _den3 * DiFfRG::za4_qcd::tr1(fenv) * _interp11 * _interp13 * _interp14 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001275510204081632, _den16 * _den3 * DiFfRG::za4_qcd::tr3(fenv) * powr<2>(_interp17) * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.002551020408163265, _den17 * _den3 * _den8 * DiFfRG::za4_qcd::tr2(fenv) * _interp11 * _interp20 * _interp21 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.01020408163265306, _den10 * _den19 * _den22 * _den4 * DiFfRG::za4_qcd::tr4(fenv) * _interp22 * _interp23 * _interp24 * _interp25 * (_interp2 * nthk3 + _interp1 * (nthk4 + 50. * (-nthk3 + nthk5))), 0.))))))))))))))))));
        // clang-format on
      }
      return _acc;
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
using DiFfRG::ZA4_kernel;