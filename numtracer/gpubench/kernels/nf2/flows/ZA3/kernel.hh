#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZA3_kernel
  {
    public:
    using Regulator = REG;

    #if NT_ZA3_QCD_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za3_qcd::nenv) > 0 ? (DiFfRG::za3_qcd::nenv) : 1];
      DiFfRG::za3_qcd::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp1 = ntRe(DiFfRG::za3_qcd::tr0(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix2 = ZA.index(l1);
      const auto _interp4 = ZA.at(_ix2);
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix5 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA.at(_ix5);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _ix4 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp8 = ZA.at(_ix4);
      const auto _ix1 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3.at(_ix1);
      const auto _ix0 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3.at(_ix0);
      const auto _ix3 = ZA3.index(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp11 = ZA3.at(_ix3);
      const auto _interp12 = ntRe(DiFfRG::za3_qcd::tr1(fenv));
      const auto _interp13 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp14 = ntRe(DiFfRG::za3_qcd::tr2(fenv));
      const auto _interp15 = ZAcbc.at(_ix1);
      const auto _interp16 = ZAcbc.at(_ix0);
      const auto _interp17 = ZAcbc.at(_ix3);
      const auto _interp18 = Zc.at(_ix2);
      const auto _interp19 = Zc.at(_ix5);
      const auto _interp20 = Zc.at(_ix4);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp4 * powr<2>(l1) + _interp2 * nthk0);
      const auto _den3 = powr<-2>(_interp18 * powr<2>(l1) + _interp2 * nthk3);
      const auto _den4 = powr<-1>(_interp5 * nthk0 + _interp6 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = -powr<-1>(_interp5 * nthk3 + _interp19 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den6 = powr<-1>(_interp7 * nthk0 + _interp8 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _den7 = -powr<-1>(_interp7 * nthk3 + _interp20 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _cse1 = powr<-2>(p);
      const auto _cse2 = _interp3 * nthk0;
      const auto _cse3 = -nthk0;
      const auto _cse4 = _cse3 + nthk2;
      const auto _cse5 = 50. * _cse4 * _den1 * powr<6>(k);
      const auto _cse6 = _cse5 + nthk1;
      const auto _cse7 = _cse6 * _interp2;
      const auto _cse8 = _cse2 + _cse7;
      return fma(-0.005050505050505051, _cse1 * _cse8 * _den2 * _den6 * _interp10 * _interp12 * _interp13, fma(-0.005050505050505051, _cse1 * _cse8 * _den2 * _den4 * _den6 * _interp1 * _interp10 * _interp11 * _interp9, fma(-0.0101010101010101, _cse1 * _den3 * _den5 * _den7 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp3 * nthk3 + _interp2 * (nthk4 + 50. * (-nthk3 + nthk5))), 0.)));
    }
    #elif NT_ZA3_QCD_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za3_qcd::nenv) > 0 ? (DiFfRG::za3_qcd::nenv) : 1];
      DiFfRG::za3_qcd::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp22 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp23 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _ix1 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp25 = ZAqbq1.at(_ix1);
      const auto _ix6 = zq.index(l1);
      const auto _interp26 = zq.at(_ix6);
      const auto _interp27 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp28 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp32 = ZAqbq4.at(_ix1);
      const auto _interp35 = ZAqbq7.at(_ix1);
      const auto _ix2 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp38 = ZAqbq1.at(_ix2);
      const auto _ix3 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp39 = ZAqbq1.at(_ix3);
      const auto _interp40 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp41 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp43 = ZAqbq4.at(_ix2);
      const auto _interp45 = ZAqbq4.at(_ix3);
      const auto _interp48 = ZAqbq7.at(_ix2);
      const auto _den2 = powr<-2>(fma(_interp22, powr<-1>(l1), _interp26));
      const auto _den5 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den6 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den8 = powr<-1>(fma(_interp40, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp41));
      const auto _den11 = powr<-1>(fma(_interp27, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))), _interp28));
      const auto _interp21 = ntIm(DiFfRG::za3_qcd::tr3(fenv));
      const auto _ix0 = ZAAqbq2.index(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp24 = ZAAqbq2.at(_ix0);
      const auto _interp29 = ntIm(DiFfRG::za3_qcd::tr4(fenv));
      const auto _interp30 = ntIm(DiFfRG::za3_qcd::tr5(fenv));
      const auto _interp31 = ZAAqbq1.at(_ix0);
      const auto _interp33 = ntIm(DiFfRG::za3_qcd::tr6(fenv));
      const auto _interp34 = ntIm(DiFfRG::za3_qcd::tr7(fenv));
      const auto _interp36 = ntIm(DiFfRG::za3_qcd::tr8(fenv));
      const auto _interp37 = ntIm(DiFfRG::za3_qcd::tr9(fenv));
      const auto _interp42 = ntIm(DiFfRG::za3_qcd::tr11(fenv));
      const auto _interp44 = ntIm(DiFfRG::za3_qcd::tr10(fenv));
      const auto _interp46 = ntIm(DiFfRG::za3_qcd::tr12(fenv));
      const auto _interp47 = ntIm(DiFfRG::za3_qcd::tr13(fenv));
      const auto _interp49 = ntIm(DiFfRG::za3_qcd::tr14(fenv));
      const auto _interp1 = ntRe(DiFfRG::za3_qcd::tr0(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp4 = ZA.at(_ix6);
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix9 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA.at(_ix9);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _ix8 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp8 = ZA.at(_ix8);
      const auto _ix5 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3.at(_ix5);
      const auto _ix4 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3.at(_ix4);
      const auto _ix7 = ZA3.index(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp11 = ZA3.at(_ix7);
      const auto _interp12 = ntRe(DiFfRG::za3_qcd::tr1(fenv));
      const auto _interp13 = ZA4.at(_ix0);
      const auto _interp14 = ntRe(DiFfRG::za3_qcd::tr2(fenv));
      const auto _interp15 = ZAcbc.at(_ix5);
      const auto _interp16 = ZAcbc.at(_ix4);
      const auto _interp17 = ZAcbc.at(_ix7);
      const auto _interp18 = Zc.at(_ix6);
      const auto _interp19 = Zc.at(_ix9);
      const auto _interp20 = Zc.at(_ix8);
      const auto _interp50 = ntIm(DiFfRG::za3_qcd::tr18(fenv));
      const auto _interp51 = ntIm(DiFfRG::za3_qcd::tr19(fenv));
      const auto _interp52 = ntIm(DiFfRG::za3_qcd::tr20(fenv));
      const auto _interp53 = ntIm(DiFfRG::za3_qcd::tr15(fenv));
      const auto _interp54 = ZAqbq7.at(_ix7);
      const auto _interp55 = ntIm(DiFfRG::za3_qcd::tr16(fenv));
      const auto _interp56 = ntIm(DiFfRG::za3_qcd::tr17(fenv));
      const auto _interp57 = ntIm(DiFfRG::za3_qcd::tr21(fenv));
      const auto _interp58 = ntIm(DiFfRG::za3_qcd::tr22(fenv));
      // clang-format off
      using _T = decltype(_den11 + _den2 + _den5 + _den6 + _den8 + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + nthk3 + nthk4 + nthk5 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _cse1_k1 = powr<-5>(l1);
        const auto _cse2_k1 = powr<-2>(p);
        const auto _cse3_k1 = -_interp22 * etaQ;
        const auto _cse4_k1 = _cse3_k1 + _interp23;
        // clang-format off
        _acc += 0.005050505050505051 * fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * _interp21 * _interp24 * _interp25, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * _interp24 * _interp25 * _interp29, fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * _interp30 * _interp31 * _interp32, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * _interp31 * _interp32 * _interp33, fma(_cse1_k1, _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * _interp24 * _interp34 * _interp35, fma(-1., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * _interp24 * _interp35 * _interp36, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp25 * _interp37 * _interp38 * _interp39, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp32 * _interp39 * _interp42 * _interp43, fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp25 * _interp43 * _interp44 * _interp45, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp32 * _interp38 * _interp45 * _interp46, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp25 * _interp39 * _interp47 * _interp48, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp32 * _interp45 * _interp48 * _interp49, 0.))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den3 = powr<-2>(_interp4 * powr<2>(l1) + _interp2 * nthk0);
        const auto _den4 = powr<-2>(_interp18 * powr<2>(l1) + _interp2 * nthk3);
        const auto _den7 = powr<-1>(_interp5 * nthk0 + _interp6 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den9 = -powr<-1>(_interp5 * nthk3 + _interp19 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp7 * nthk0 + _interp8 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp7 * nthk3 + _interp20 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1_k2 = powr<-2>(p);
        const auto _cse2_k2 = _interp3 * nthk0;
        const auto _cse3_k2 = -nthk0;
        const auto _cse4_k2 = _cse3_k2 + nthk2;
        const auto _cse5_k2 = 50. * _cse4_k2 * _den1 * powr<6>(k);
        const auto _cse6_k2 = _cse5_k2 + nthk1;
        const auto _cse7_k2 = _cse6_k2 * _interp2;
        const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
        // clang-format off
        _acc += fma(-0.005050505050505051, _cse1_k2 * _cse8_k2 * _den10 * _den3 * _interp10 * _interp12 * _interp13, fma(-0.005050505050505051, _cse1_k2 * _cse8_k2 * _den10 * _den3 * _den7 * _interp1 * _interp10 * _interp11 * _interp9, fma(-0.0101010101010101, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp35 * _interp38 * _interp39 * _interp50 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(0.0101010101010101, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp35 * _interp43 * _interp45 * _interp51 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.005050505050505051, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp35 * _interp39 * _interp48 * _interp52 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.0101010101010101, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp25 * _interp38 * _interp53 * _interp54 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.0101010101010101, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp32 * _interp43 * _interp54 * _interp55 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.005050505050505051, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp25 * _interp48 * _interp54 * _interp56 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.005050505050505051, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp35 * _interp38 * _interp54 * _interp57 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.002525252525252526, _cse1_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * _interp35 * _interp48 * _interp54 * _interp58 * (_interp23 - _interp22 * etaQ) * powr<-5>(l1), fma(-0.0101010101010101, _cse1_k2 * _den12 * _den4 * _den9 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp3 * nthk3 + _interp2 * (nthk4 + 50. * (-nthk3 + nthk5))), 0.)))))))))));
        // clang-format on
      }
      return _acc;
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za3_qcd::nenv) > 0 ? (DiFfRG::za3_qcd::nenv) : 1];
      DiFfRG::za3_qcd::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp18 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp19 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _ix1 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp21 = ZAqbq1.at(_ix1);
      const auto _ix6 = zq.index(l1);
      const auto _interp22 = zq.at(_ix6);
      const auto _interp23 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp24 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp26 = ZAqbq4.at(_ix1);
      const auto _interp27 = ZAqbq7.at(_ix1);
      const auto _ix2 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp28 = ZAqbq1.at(_ix2);
      const auto _ix3 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp29 = ZAqbq1.at(_ix3);
      const auto _interp30 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp31 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = ZAqbq4.at(_ix2);
      const auto _interp33 = ZAqbq4.at(_ix3);
      const auto _den2 = powr<-2>(fma(_interp18, powr<-1>(l1), _interp22));
      const auto _den5 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den6 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den8 = powr<-1>(fma(_interp30, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp31));
      const auto _den11 = powr<-1>(fma(_interp23, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))), _interp24));
      const auto _ix0 = ZAAqbq2.index(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp20 = ZAAqbq2.at(_ix0);
      const auto _interp25 = ZAAqbq1.at(_ix0);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA.at(_ix6);
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix9 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp5 = ZA.at(_ix9);
      const auto _interp6 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _ix8 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp7 = ZA.at(_ix8);
      const auto _ix5 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp8 = ZA3.at(_ix5);
      const auto _ix4 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3.at(_ix4);
      const auto _ix7 = ZA3.index(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3.at(_ix7);
      const auto _interp11 = ZA4.at(_ix0);
      const auto _interp12 = ZAcbc.at(_ix5);
      const auto _interp13 = ZAcbc.at(_ix4);
      const auto _interp14 = ZAcbc.at(_ix7);
      const auto _interp15 = Zc.at(_ix6);
      const auto _interp16 = Zc.at(_ix9);
      const auto _interp17 = Zc.at(_ix8);
      const auto _interp34 = ZAqbq7.at(_ix5);
      const auto _interp35 = ZAqbq7.at(_ix7);
      // clang-format off
      using _T = decltype(complex<double>(0.,0.0101010101010101) + _den11 + _den2 + _den5 + _den6 + _den8 + DiFfRG::za3_qcd::tr0(fenv) + DiFfRG::za3_qcd::tr10(fenv) + DiFfRG::za3_qcd::tr11(fenv) + DiFfRG::za3_qcd::tr12(fenv) + DiFfRG::za3_qcd::tr13(fenv) + DiFfRG::za3_qcd::tr14(fenv) + DiFfRG::za3_qcd::tr15(fenv) + DiFfRG::za3_qcd::tr16(fenv) + DiFfRG::za3_qcd::tr17(fenv) + DiFfRG::za3_qcd::tr18(fenv) + DiFfRG::za3_qcd::tr19(fenv) + DiFfRG::za3_qcd::tr1(fenv) + DiFfRG::za3_qcd::tr20(fenv) + DiFfRG::za3_qcd::tr21(fenv) + DiFfRG::za3_qcd::tr22(fenv) + DiFfRG::za3_qcd::tr2(fenv) + DiFfRG::za3_qcd::tr3(fenv) + DiFfRG::za3_qcd::tr4(fenv) + DiFfRG::za3_qcd::tr5(fenv) + DiFfRG::za3_qcd::tr6(fenv) + DiFfRG::za3_qcd::tr7(fenv) + DiFfRG::za3_qcd::tr8(fenv) + DiFfRG::za3_qcd::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp4 + _interp5 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + nthk3 + nthk4 + nthk5 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _cse1_k1 = powr<-5>(l1);
        const auto _cse2_k1 = powr<-2>(p);
        const auto _cse3_k1 = -_interp18 * etaQ;
        const auto _cse4_k1 = _cse3_k1 + _interp19;
        // clang-format off
        _acc += complex<double>(0.,-0.005050505050505051) * fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * DiFfRG::za3_qcd::tr3(fenv) * _interp20 * _interp21, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * DiFfRG::za3_qcd::tr4(fenv) * _interp20 * _interp21, fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * DiFfRG::za3_qcd::tr5(fenv) * _interp25 * _interp26, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * DiFfRG::za3_qcd::tr6(fenv) * _interp25 * _interp26, fma(_cse1_k1, _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * DiFfRG::za3_qcd::tr7(fenv) * _interp20 * _interp27, fma(-1., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den6 * DiFfRG::za3_qcd::tr8(fenv) * _interp20 * _interp27, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr9(fenv) * _interp21 * _interp28 * _interp29, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr18(fenv) * _interp27 * _interp28 * _interp29, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr11(fenv) * _interp26 * _interp29 * _interp32, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr12(fenv) * _interp26 * _interp28 * _interp33, fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr10(fenv) * _interp21 * _interp32 * _interp33, fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr19(fenv) * _interp27 * _interp32 * _interp33, 0.))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den3 = powr<-2>(_interp3 * powr<2>(l1) + _interp1 * nthk0);
        const auto _den4 = powr<-2>(_interp15 * powr<2>(l1) + _interp1 * nthk3);
        const auto _den7 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den9 = -powr<-1>(_interp4 * nthk3 + _interp16 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp6 * nthk0 + _interp7 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp6 * nthk3 + _interp17 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1_k2 = powr<-2>(p);
        const auto _cse2_k2 = _interp2 * nthk0;
        const auto _cse3_k2 = -nthk0;
        const auto _cse4_k2 = _cse3_k2 + nthk2;
        const auto _cse5_k2 = 50. * _cse4_k2 * _den1 * powr<6>(k);
        const auto _cse6_k2 = _cse5_k2 + nthk1;
        const auto _cse7_k2 = _cse6_k2 * _interp1;
        const auto _cse8_k2 = _cse2_k2 + _cse7_k2;
        const auto _cse9_k2 = powr<-5>(l1);
        const auto _cse10_k2 = -_interp18 * etaQ;
        const auto _cse11_k2 = _cse10_k2 + _interp19;
        // clang-format off
        _acc += fma(complex<double>(0.,0.0101010101010101), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr13(fenv) * _interp21 * _interp29 * _interp34, fma(complex<double>(0.,0.005050505050505051), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr20(fenv) * _interp27 * _interp29 * _interp34, fma(complex<double>(0.,0.0101010101010101), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr14(fenv) * _interp26 * _interp33 * _interp34, fma(complex<double>(0.,0.0101010101010101), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr15(fenv) * _interp21 * _interp28 * _interp35, fma(complex<double>(0.,0.005050505050505051), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr21(fenv) * _interp27 * _interp28 * _interp35, fma(complex<double>(0.,0.0101010101010101), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr16(fenv) * _interp26 * _interp32 * _interp35, fma(complex<double>(0.,0.005050505050505051), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr17(fenv) * _interp21 * _interp34 * _interp35, fma(complex<double>(0.,0.002525252525252526), _cse11_k2 * _cse1_k2 * _cse9_k2 * _den11 * _den2 * _den5 * _den6 * _den8 * DiFfRG::za3_qcd::tr22(fenv) * _interp27 * _interp34 * _interp35, fma(-0.005050505050505051, _cse1_k2 * _cse8_k2 * _den10 * _den3 * DiFfRG::za3_qcd::tr1(fenv) * _interp11 * _interp9, fma(-0.005050505050505051, _cse1_k2 * _cse8_k2 * _den10 * _den3 * _den7 * DiFfRG::za3_qcd::tr0(fenv) * _interp10 * _interp8 * _interp9, fma(-0.0101010101010101, _cse1_k2 * _den12 * _den4 * _den9 * DiFfRG::za3_qcd::tr2(fenv) * _interp12 * _interp13 * _interp14 * (_interp2 * nthk3 + _interp1 * (nthk4 + 50. * (-nthk3 + nthk5))), 0.)))))))))));
        // clang-format on
      }
      return _acc;
    }
    #endif

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

    static KOKKOS_INLINE_FUNCTION double ntRe(double x) { return x; }
    template <class T> static KOKKOS_INLINE_FUNCTION double ntRe(const T &z) { return z.real(); }
    static KOKKOS_INLINE_FUNCTION double ntIm(double) { return 0.0; }
    template <class T> static KOKKOS_INLINE_FUNCTION double ntIm(const T &z) { return z.imag(); }
  };
}
using DiFfRG::ZA3_kernel;