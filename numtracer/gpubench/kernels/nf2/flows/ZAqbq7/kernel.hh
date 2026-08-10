#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZAqbq7_kernel
  {
    public:
    using Regulator = REG;

    #if NT_ZAQBQ7_QCD_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaqbq7_qcd::nenv) > 0 ? (DiFfRG::zaqbq7_qcd::nenv) : 1];
      DiFfRG::zaqbq7_qcd::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp2 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp4 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _ix11 = ZA.index(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp5 = ZA.at(_ix11);
      const auto _ix1 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1)))));
      const auto _interp6 = ZAAqbq2.at(_ix1);
      const auto _ix5 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp7 = ZAqbq1.at(_ix5);
      const auto _ix8 = zq.index(l1);
      const auto _interp8 = zq.at(_ix8);
      const auto _interp10 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _ix10 = ZA.index(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp11 = ZA.at(_ix10);
      const auto _ix0 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp12 = ZAAqbq2.at(_ix0);
      const auto _ix2 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp13 = ZAqbq1.at(_ix2);
      const auto _interp15 = ZAAqbq1.at(_ix1);
      const auto _interp16 = ZAqbq4.at(_ix5);
      const auto _interp18 = ZAAqbq1.at(_ix0);
      const auto _interp19 = ZAqbq4.at(_ix2);
      const auto _interp21 = ZAqbq7.at(_ix5);
      const auto _interp23 = ZAqbq7.at(_ix2);
      const auto _ix4 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp25 = ZAqbq1.at(_ix4);
      const auto _ix6 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp26 = ZAqbq1.at(_ix6);
      const auto _interp27 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp28 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp30 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp31 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _ix7 = ZAqbq1.index(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = ZAqbq1.at(_ix7);
      const auto _ix3 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp33 = ZAqbq1.at(_ix3);
      const auto _interp35 = ZAqbq4.at(_ix4);
      const auto _interp37 = ZAqbq4.at(_ix7);
      const auto _interp39 = ZAqbq4.at(_ix6);
      const auto _interp42 = ZAqbq4.at(_ix3);
      const auto _interp45 = ZAqbq7.at(_ix4);
      const auto _interp53 = ZAqbq7.at(_ix7);
      const auto _interp63 = ZAqbq7.at(_ix3);
      const auto _interp69 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp70 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp71 = ZA.at(_ix8);
      const auto _interp72 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp73 = zq.at(_ix11);
      const auto _interp77 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp78 = zq.at(_ix10);
      const auto _interp82 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp2, powr<-1>(l1), _interp8));
      const auto _den3 = powr<-2>(fma(_interp71, powr<2>(l1), fma(_interp69, nthk0, 0.)));
      const auto _den4 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den5 = powr<-1>(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den6 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den7 = powr<-1>(fma(_interp4, nthk0, fma(_interp5, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(_interp30, nthk0, fma(_interp31, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den9 = powr<-1>(fma(_interp27, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp28));
      const auto _den10 = powr<-1>(fma(_interp72, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))), _interp73));
      const auto _den11 = powr<-1>(fma(_interp10, nthk0, fma(_interp11, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den12 = powr<-1>(fma(_interp77, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))), _interp78));
      const auto _interp14 = ntRe(DiFfRG::zaqbq7_qcd::tr4(fenv));
      const auto _interp17 = ntRe(DiFfRG::zaqbq7_qcd::tr5(fenv));
      const auto _interp22 = ntRe(DiFfRG::zaqbq7_qcd::tr9(fenv));
      const auto _interp24 = ntRe(DiFfRG::zaqbq7_qcd::tr13(fenv));
      const auto _interp29 = ntRe(DiFfRG::zaqbq7_qcd::tr12(fenv));
      const auto _interp34 = ntRe(DiFfRG::zaqbq7_qcd::tr18(fenv));
      const auto _interp36 = ntRe(DiFfRG::zaqbq7_qcd::tr14(fenv));
      const auto _interp38 = ntRe(DiFfRG::zaqbq7_qcd::tr15(fenv));
      const auto _interp40 = ntRe(DiFfRG::zaqbq7_qcd::tr19(fenv));
      const auto _interp41 = ntRe(DiFfRG::zaqbq7_qcd::tr16(fenv));
      const auto _interp43 = ntRe(DiFfRG::zaqbq7_qcd::tr17(fenv));
      const auto _interp44 = ntRe(DiFfRG::zaqbq7_qcd::tr21(fenv));
      const auto _interp46 = ntRe(DiFfRG::zaqbq7_qcd::tr20(fenv));
      const auto _interp47 = ntRe(DiFfRG::zaqbq7_qcd::tr23(fenv));
      const auto _interp20 = ntRe(DiFfRG::zaqbq7_qcd::tr7(fenv));
      const auto _interp48 = ntRe(DiFfRG::zaqbq7_qcd::tr22(fenv));
      const auto _interp49 = ntRe(DiFfRG::zaqbq7_qcd::tr35(fenv));
      const auto _interp50 = ntRe(DiFfRG::zaqbq7_qcd::tr36(fenv));
      const auto _interp51 = ntRe(DiFfRG::zaqbq7_qcd::tr37(fenv));
      const auto _interp52 = ntRe(DiFfRG::zaqbq7_qcd::tr24(fenv));
      const auto _interp54 = ntRe(DiFfRG::zaqbq7_qcd::tr25(fenv));
      const auto _interp55 = ntRe(DiFfRG::zaqbq7_qcd::tr26(fenv));
      const auto _interp56 = ntRe(DiFfRG::zaqbq7_qcd::tr27(fenv));
      const auto _interp57 = ZAqbq7(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp58 = ntRe(DiFfRG::zaqbq7_qcd::tr28(fenv));
      const auto _interp59 = ntRe(DiFfRG::zaqbq7_qcd::tr29(fenv));
      const auto _interp60 = ntRe(DiFfRG::zaqbq7_qcd::tr38(fenv));
      const auto _interp61 = ntRe(DiFfRG::zaqbq7_qcd::tr39(fenv));
      const auto _interp62 = ntRe(DiFfRG::zaqbq7_qcd::tr30(fenv));
      const auto _interp1 = ntRe(DiFfRG::zaqbq7_qcd::tr1(fenv));
      const auto _interp9 = ntRe(DiFfRG::zaqbq7_qcd::tr2(fenv));
      const auto _interp64 = ntRe(DiFfRG::zaqbq7_qcd::tr31(fenv));
      const auto _interp65 = ntRe(DiFfRG::zaqbq7_qcd::tr32(fenv));
      const auto _interp66 = ntRe(DiFfRG::zaqbq7_qcd::tr33(fenv));
      const auto _interp67 = ntRe(DiFfRG::zaqbq7_qcd::tr34(fenv));
      const auto _interp68 = ntRe(DiFfRG::zaqbq7_qcd::tr40(fenv));
      const auto _interp74 = ntRe(DiFfRG::zaqbq7_qcd::tr41(fenv));
      const auto _interp75 = ntRe(DiFfRG::zaqbq7_qcd::tr42(fenv));
      const auto _interp76 = ntRe(DiFfRG::zaqbq7_qcd::tr48(fenv));
      const auto _interp79 = ntRe(DiFfRG::zaqbq7_qcd::tr50(fenv));
      const auto _interp80 = ntRe(DiFfRG::zaqbq7_qcd::tr53(fenv));
      const auto _interp81 = ntRe(DiFfRG::zaqbq7_qcd::tr56(fenv));
      const auto _interp83 = ntRe(DiFfRG::zaqbq7_qcd::tr58(fenv));
      const auto _interp84 = ntRe(DiFfRG::zaqbq7_qcd::tr57(fenv));
      const auto _ix9 = ZAqbq4.index(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp85 = ZAqbq4.at(_ix9);
      const auto _interp86 = ntRe(DiFfRG::zaqbq7_qcd::tr59(fenv));
      const auto _interp87 = ntRe(DiFfRG::zaqbq7_qcd::tr60(fenv));
      const auto _interp88 = ntRe(DiFfRG::zaqbq7_qcd::tr61(fenv));
      const auto _interp89 = ntRe(DiFfRG::zaqbq7_qcd::tr65(fenv));
      const auto _interp90 = ntRe(DiFfRG::zaqbq7_qcd::tr66(fenv));
      const auto _interp91 = ntRe(DiFfRG::zaqbq7_qcd::tr67(fenv));
      const auto _interp92 = ntRe(DiFfRG::zaqbq7_qcd::tr62(fenv));
      const auto _interp93 = ZAqbq7.at(_ix9);
      const auto _interp94 = ntRe(DiFfRG::zaqbq7_qcd::tr63(fenv));
      const auto _interp95 = ntRe(DiFfRG::zaqbq7_qcd::tr64(fenv));
      const auto _interp96 = ntRe(DiFfRG::zaqbq7_qcd::tr68(fenv));
      const auto _interp97 = ntRe(DiFfRG::zaqbq7_qcd::tr69(fenv));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den11 + _den12 + _den2 + _den3 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
      // clang-format off
      _acc += 0.01041666666666667 * fma(2., _den2 * _den7 * _interp14 * _interp15 * _interp16 * _interp3 * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _interp17 * _interp18 * _interp19 * _interp3 * powr<-5>(l1) * powr<-4>(p), fma(_den11, _den2 * _interp12 * _interp22 * _interp23 * _interp3 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp24 * _interp25 * _interp26 * _interp3 * powr<-5>(l1) * powr<-4>(p), fma(2., _den2 * _den4 * _den8 * _den9 * _interp25 * _interp29 * _interp3 * _interp32 * _interp33 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp26 * _interp3 * _interp34 * _interp35 * powr<-5>(l1) * powr<-4>(p), fma(2., _den2 * _den4 * _den8 * _den9 * _interp3 * _interp33 * _interp35 * _interp36 * _interp37 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp3 * _interp35 * _interp38 * _interp39 * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp25 * _interp3 * _interp39 * _interp40 * powr<-5>(l1) * powr<-4>(p), fma(2., _den2 * _den4 * _den8 * _den9 * _interp3 * _interp32 * _interp35 * _interp41 * _interp42 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den4 * _den8 * _den9 * _interp25 * _interp3 * _interp37 * _interp42 * _interp43 * powr<-5>(l1) * powr<-4>(p), fma(-1., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp26 * _interp3 * _interp44 * _interp45 * powr<-5>(l1) * powr<-4>(p), fma(_den2, _den4 * _den8 * _den9 * _interp3 * _interp32 * _interp33 * _interp45 * _interp46 * powr<-5>(l1) * powr<-4>(p), fma(_den11, _den2 * _den4 * _den9 * _interp19 * _interp3 * _interp39 * _interp45 * _interp47 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den7 * _interp14 * _interp15 * _interp16 * _interp2 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _interp17 * _interp18 * _interp19 * _interp2 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-1., _den11 * _den2 * _interp12 * _interp2 * _interp22 * _interp23 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp2 * _interp24 * _interp25 * _interp26 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp25 * _interp29 * _interp32 * _interp33 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp2 * _interp26 * _interp34 * _interp35 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp33 * _interp35 * _interp36 * _interp37 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp2 * _interp35 * _interp38 * _interp39 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp2 * _interp25 * _interp39 * _interp40 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp32 * _interp35 * _interp41 * _interp42 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp25 * _interp37 * _interp42 * _interp43 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(_den11, _den2 * _den4 * _den9 * _interp13 * _interp2 * _interp26 * _interp44 * _interp45 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-1., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp32 * _interp33 * _interp45 * _interp46 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-1., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp2 * _interp39 * _interp45 * _interp47 * etaQ * powr<-5>(l1) * powr<-4>(p), 0.))))))))))))))))))))))))))));
      // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += 0.002604166666666666 * fma(-4., _den2 * _den4 * _den8 * _den9 * _interp3 * _interp37 * _interp42 * _interp45 * _interp48 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den2 * _den4 * _den9 * _interp23 * _interp25 * _interp26 * _interp3 * _interp49 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den2 * _den4 * _den9 * _interp23 * _interp3 * _interp35 * _interp39 * _interp50 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp23 * _interp26 * _interp3 * _interp45 * _interp51 * powr<-5>(l1) * powr<-4>(p), fma(4., _den2 * _den4 * _den8 * _den9 * _interp25 * _interp3 * _interp33 * _interp52 * _interp53 * powr<-5>(l1) * powr<-4>(p), fma(4., _den2 * _den4 * _den8 * _den9 * _interp3 * _interp35 * _interp42 * _interp53 * _interp54 * powr<-5>(l1) * powr<-4>(p), fma(2., _den2 * _den4 * _den8 * _den9 * _interp3 * _interp33 * _interp45 * _interp53 * _interp55 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp25 * _interp3 * _interp56 * _interp57 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp3 * _interp35 * _interp57 * _interp58 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp3 * _interp45 * _interp57 * _interp59 * powr<-5>(l1) * powr<-4>(p), fma(4., _den2 * _den7 * _interp20 * _interp21 * _interp3 * _interp6 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den2 * _den4 * _den9 * _interp23 * _interp25 * _interp3 * _interp57 * _interp60 * powr<-5>(l1) * powr<-4>(p), fma(-1., _den11 * _den2 * _den4 * _den9 * _interp23 * _interp3 * _interp45 * _interp57 * _interp61 * powr<-5>(l1) * powr<-4>(p), fma(4., _den2 * _den4 * _den8 * _den9 * _interp25 * _interp3 * _interp32 * _interp62 * _interp63 * powr<-5>(l1) * powr<-4>(p), fma(4., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp37 * _interp42 * _interp45 * _interp48 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den2 * _den4 * _den9 * _interp2 * _interp23 * _interp25 * _interp26 * _interp49 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den2 * _den4 * _den9 * _interp2 * _interp23 * _interp35 * _interp39 * _interp50 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp2 * _interp23 * _interp26 * _interp45 * _interp51 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp25 * _interp33 * _interp52 * _interp53 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp35 * _interp42 * _interp53 * _interp54 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp33 * _interp45 * _interp53 * _interp55 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp2 * _interp25 * _interp56 * _interp57 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den2 * _den4 * _den9 * _interp19 * _interp2 * _interp35 * _interp57 * _interp58 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp13 * _interp2 * _interp45 * _interp57 * _interp59 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den2 * _den7 * _interp2 * _interp20 * _interp21 * _interp6 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den2 * _den4 * _den9 * _interp2 * _interp23 * _interp25 * _interp57 * _interp60 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(_den11, _den2 * _den4 * _den9 * _interp2 * _interp23 * _interp45 * _interp57 * _interp61 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den2 * _den4 * _den8 * _den9 * _interp2 * _interp25 * _interp32 * _interp62 * _interp63 * etaQ * powr<-5>(l1) * powr<-4>(p), 0.))))))))))))))))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.01041666666666667, _den2 * _den4 * _den8 * _den9 * _interp35 * _interp37 * _interp63 * _interp64 * (_interp3 - _interp2 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.005208333333333332, _den2 * _den4 * _den8 * _den9 * _interp32 * _interp45 * _interp63 * _interp65 * (_interp3 - _interp2 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.005208333333333332, _den2 * _den4 * _den8 * _den9 * _interp25 * _interp53 * _interp63 * _interp66 * (_interp3 - _interp2 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.002604166666666666, _den2 * _den4 * _den8 * _den9 * _interp45 * _interp53 * _interp63 * _interp67 * (_interp3 - _interp2 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.02083333333333333, _den2 * _den7 * _interp1 * _interp6 * _interp7 * (_interp3 - _interp2 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.02083333333333333, _den11 * _den2 * _interp12 * _interp13 * _interp9 * (_interp3 - _interp2 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.02083333333333333, _den10 * _den3 * _den5 * _interp6 * _interp68 * _interp7 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den10 * _den3 * _den5 * _interp15 * _interp16 * _interp74 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den10 * _den3 * _den5 * _interp21 * _interp6 * _interp75 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den12 * _den3 * _den6 * _interp12 * _interp13 * _interp76 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den12 * _den3 * _den6 * _interp18 * _interp19 * _interp79 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den12 * _den3 * _den6 * _interp12 * _interp23 * _interp80 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den10 * _den12 * _den3 * _den5 * _den6 * _interp13 * _interp7 * _interp81 * _interp82 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den10 * _den12 * _den3 * _den5 * _den6 * _interp16 * _interp19 * _interp82 * _interp83 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(-0.02083333333333333, _den10 * _den12 * _den3 * _den5 * _den6 * _interp13 * _interp16 * _interp84 * _interp85 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den10 * _den12 * _den3 * _den5 * _den6 * _interp19 * _interp7 * _interp85 * _interp86 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den10 * _den12 * _den3 * _den5 * _den6 * _interp13 * _interp21 * _interp82 * _interp87 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den10 * _den12 * _den3 * _den5 * _den6 * _interp19 * _interp21 * _interp85 * _interp88 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den10 * _den12 * _den3 * _den5 * _den6 * _interp23 * _interp7 * _interp82 * _interp89 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den10 * _den12 * _den3 * _den5 * _den6 * _interp16 * _interp23 * _interp85 * _interp90 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.005208333333333332, _den10 * _den12 * _den3 * _den5 * _den6 * _interp21 * _interp23 * _interp82 * _interp91 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den10 * _den12 * _den3 * _den5 * _den6 * _interp13 * _interp7 * _interp92 * _interp93 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den10 * _den12 * _den3 * _den5 * _den6 * _interp16 * _interp19 * _interp93 * _interp94 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.005208333333333332, _den10 * _den12 * _den3 * _den5 * _den6 * _interp13 * _interp21 * _interp93 * _interp95 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.005208333333333332, _den10 * _den12 * _den3 * _den5 * _den6 * _interp23 * _interp7 * _interp93 * _interp96 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.002604166666666666, _den10 * _den12 * _den3 * _den5 * _den6 * _interp21 * _interp23 * _interp93 * _interp97 * (_interp70 * nthk0 + _interp69 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), 0.))))))))))));
      // clang-format on
      }
      return _acc;
    }
    #elif NT_ZAQBQ7_QCD_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaqbq7_qcd::nenv) > 0 ? (DiFfRG::zaqbq7_qcd::nenv) : 1];
      DiFfRG::zaqbq7_qcd::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix9 = ZA.index(l1);
      const auto _interp4 = ZA.at(_ix9);
      const auto _interp5 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _ix12 = ZA.index(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp6 = ZA.at(_ix12);
      const auto _ix5 = ZA3.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp7 = ZA3.at(_ix5);
      const auto _interp10 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp11 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp12 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _ix13 = ZA.index(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp13 = ZA.at(_ix13);
      const auto _ix1 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1)))));
      const auto _interp14 = ZAAqbq2.at(_ix1);
      const auto _ix6 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp15 = ZAqbq1.at(_ix6);
      const auto _interp16 = zq.at(_ix9);
      const auto _interp18 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _ix11 = ZA.index(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp19 = ZA.at(_ix11);
      const auto _ix0 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp20 = ZAAqbq2.at(_ix0);
      const auto _ix2 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp21 = ZAqbq1.at(_ix2);
      const auto _interp25 = ZAAqbq1.at(_ix1);
      const auto _interp26 = ZAqbq4.at(_ix6);
      const auto _interp28 = ZAAqbq1.at(_ix0);
      const auto _interp29 = ZAqbq4.at(_ix2);
      const auto _interp32 = ZAqbq7.at(_ix6);
      const auto _interp35 = ZAqbq7.at(_ix2);
      const auto _interp39 = ZAqbq1.at(_ix5);
      const auto _ix7 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp40 = ZAqbq1.at(_ix7);
      const auto _interp41 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp42 = zq.at(_ix12);
      const auto _interp44 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp45 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _ix8 = ZAqbq1.index(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp46 = ZAqbq1.at(_ix8);
      const auto _ix4 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp47 = ZAqbq1.at(_ix4);
      const auto _interp49 = ZAqbq4.at(_ix5);
      const auto _interp51 = ZAqbq4.at(_ix8);
      const auto _interp53 = ZAqbq4.at(_ix7);
      const auto _interp56 = ZAqbq4.at(_ix4);
      const auto _interp59 = ZAqbq7.at(_ix5);
      const auto _interp67 = ZAqbq7.at(_ix8);
      const auto _interp71 = ZAqbq7.at(_ix7);
      const auto _interp77 = ZAqbq7.at(_ix4);
      const auto _interp83 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp84 = zq.at(_ix13);
      const auto _interp88 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp89 = zq.at(_ix11);
      const auto _ix3 = ZAqbq4.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp101 = ZAqbq4.at(_ix3);
      const auto _interp109 = ZAqbq7.at(_ix3);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp10, powr<-1>(l1), _interp16));
      const auto _den3 = powr<-2>(fma(_interp4, powr<2>(l1), fma(_interp2, nthk0, 0.)));
      const auto _den4 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den5 = powr<-1>(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den7 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den8 = powr<-1>(fma(_interp5, nthk0, fma(_interp6, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den9 = powr<-1>(fma(_interp12, nthk0, fma(_interp13, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den10 = powr<-1>(fma(_interp44, nthk0, fma(_interp45, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den11 = powr<-1>(fma(_interp41, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp42));
      const auto _den12 = powr<-1>(fma(_interp83, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))), _interp84));
      const auto _den14 = powr<-1>(fma(_interp18, nthk0, fma(_interp19, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den15 = powr<-1>(fma(_interp88, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))), _interp89));
      const auto _interp17 = ntRe(DiFfRG::zaqbq7_qcd::tr2(fenv));
      const auto _interp22 = ntIm(DiFfRG::zaqbq7_qcd::tr3(fenv));
      const auto _ix10 = ZA3.index(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp23 = ZA3.at(_ix10);
      const auto _interp24 = ntRe(DiFfRG::zaqbq7_qcd::tr4(fenv));
      const auto _interp27 = ntRe(DiFfRG::zaqbq7_qcd::tr5(fenv));
      const auto _interp30 = ntIm(DiFfRG::zaqbq7_qcd::tr6(fenv));
      const auto _interp31 = ntRe(DiFfRG::zaqbq7_qcd::tr7(fenv));
      const auto _interp33 = ntIm(DiFfRG::zaqbq7_qcd::tr8(fenv));
      const auto _interp34 = ntRe(DiFfRG::zaqbq7_qcd::tr9(fenv));
      const auto _interp36 = ntIm(DiFfRG::zaqbq7_qcd::tr10(fenv));
      const auto _interp37 = ntIm(DiFfRG::zaqbq7_qcd::tr11(fenv));
      const auto _interp38 = ntRe(DiFfRG::zaqbq7_qcd::tr13(fenv));
      const auto _interp43 = ntRe(DiFfRG::zaqbq7_qcd::tr12(fenv));
      const auto _interp48 = ntRe(DiFfRG::zaqbq7_qcd::tr18(fenv));
      const auto _interp50 = ntRe(DiFfRG::zaqbq7_qcd::tr14(fenv));
      const auto _interp52 = ntRe(DiFfRG::zaqbq7_qcd::tr15(fenv));
      const auto _interp54 = ntRe(DiFfRG::zaqbq7_qcd::tr19(fenv));
      const auto _interp55 = ntRe(DiFfRG::zaqbq7_qcd::tr16(fenv));
      const auto _interp57 = ntRe(DiFfRG::zaqbq7_qcd::tr17(fenv));
      const auto _interp58 = ntRe(DiFfRG::zaqbq7_qcd::tr21(fenv));
      const auto _interp60 = ntRe(DiFfRG::zaqbq7_qcd::tr20(fenv));
      const auto _interp61 = ntRe(DiFfRG::zaqbq7_qcd::tr23(fenv));
      const auto _interp62 = ntRe(DiFfRG::zaqbq7_qcd::tr22(fenv));
      const auto _interp63 = ntRe(DiFfRG::zaqbq7_qcd::tr35(fenv));
      const auto _interp64 = ntRe(DiFfRG::zaqbq7_qcd::tr36(fenv));
      const auto _interp65 = ntRe(DiFfRG::zaqbq7_qcd::tr37(fenv));
      const auto _interp66 = ntRe(DiFfRG::zaqbq7_qcd::tr24(fenv));
      const auto _interp68 = ntRe(DiFfRG::zaqbq7_qcd::tr25(fenv));
      const auto _interp69 = ntRe(DiFfRG::zaqbq7_qcd::tr26(fenv));
      const auto _interp9 = ntRe(DiFfRG::zaqbq7_qcd::tr1(fenv));
      const auto _interp70 = ntRe(DiFfRG::zaqbq7_qcd::tr27(fenv));
      const auto _interp72 = ntRe(DiFfRG::zaqbq7_qcd::tr28(fenv));
      const auto _interp73 = ntRe(DiFfRG::zaqbq7_qcd::tr29(fenv));
      const auto _interp74 = ntRe(DiFfRG::zaqbq7_qcd::tr38(fenv));
      const auto _interp75 = ntRe(DiFfRG::zaqbq7_qcd::tr39(fenv));
      const auto _interp76 = ntRe(DiFfRG::zaqbq7_qcd::tr30(fenv));
      const auto _interp78 = ntRe(DiFfRG::zaqbq7_qcd::tr31(fenv));
      const auto _interp79 = ntRe(DiFfRG::zaqbq7_qcd::tr32(fenv));
      const auto _interp80 = ntRe(DiFfRG::zaqbq7_qcd::tr33(fenv));
      const auto _interp81 = ntRe(DiFfRG::zaqbq7_qcd::tr34(fenv));
      const auto _interp100 = ntRe(DiFfRG::zaqbq7_qcd::tr57(fenv));
      const auto _interp102 = ntRe(DiFfRG::zaqbq7_qcd::tr59(fenv));
      const auto _interp108 = ntRe(DiFfRG::zaqbq7_qcd::tr62(fenv));
      const auto _interp1 = ntIm(DiFfRG::zaqbq7_qcd::tr0(fenv));
      const auto _interp8 = ZAAqbq2(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
      const auto _interp82 = ntRe(DiFfRG::zaqbq7_qcd::tr40(fenv));
      const auto _interp85 = ntRe(DiFfRG::zaqbq7_qcd::tr41(fenv));
      const auto _interp104 = ntRe(DiFfRG::zaqbq7_qcd::tr61(fenv));
      const auto _interp106 = ntRe(DiFfRG::zaqbq7_qcd::tr66(fenv));
      const auto _interp110 = ntRe(DiFfRG::zaqbq7_qcd::tr63(fenv));
      const auto _interp111 = ntRe(DiFfRG::zaqbq7_qcd::tr64(fenv));
      const auto _interp112 = ntRe(DiFfRG::zaqbq7_qcd::tr68(fenv));
      const auto _interp113 = ntRe(DiFfRG::zaqbq7_qcd::tr69(fenv));
      const auto _interp114 = ntIm(DiFfRG::zaqbq7_qcd::tr43(fenv));
      const auto _interp115 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp116 = zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp117 = ntIm(DiFfRG::zaqbq7_qcd::tr44(fenv));
      const auto _interp118 = ntIm(DiFfRG::zaqbq7_qcd::tr45(fenv));
      const auto _interp119 = ntIm(DiFfRG::zaqbq7_qcd::tr46(fenv));
      const auto _interp120 = ntIm(DiFfRG::zaqbq7_qcd::tr47(fenv));
      const auto _den13 = powr<-1>(_interp116 + _interp115 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _interp86 = ntRe(DiFfRG::zaqbq7_qcd::tr42(fenv));
      const auto _interp87 = ntRe(DiFfRG::zaqbq7_qcd::tr48(fenv));
      const auto _interp90 = ntIm(DiFfRG::zaqbq7_qcd::tr49(fenv));
      const auto _interp91 = ntRe(DiFfRG::zaqbq7_qcd::tr50(fenv));
      const auto _interp92 = ntIm(DiFfRG::zaqbq7_qcd::tr51(fenv));
      const auto _interp93 = ntRe(DiFfRG::zaqbq7_qcd::tr53(fenv));
      const auto _interp94 = ntIm(DiFfRG::zaqbq7_qcd::tr54(fenv));
      const auto _interp95 = ntIm(DiFfRG::zaqbq7_qcd::tr52(fenv));
      const auto _interp96 = ntIm(DiFfRG::zaqbq7_qcd::tr55(fenv));
      const auto _interp97 = ntRe(DiFfRG::zaqbq7_qcd::tr56(fenv));
      const auto _interp98 = ZAqbq1.at(_ix10);
      const auto _interp99 = ntRe(DiFfRG::zaqbq7_qcd::tr58(fenv));
      const auto _interp103 = ntRe(DiFfRG::zaqbq7_qcd::tr60(fenv));
      const auto _interp105 = ntRe(DiFfRG::zaqbq7_qcd::tr65(fenv));
      const auto _interp107 = ntRe(DiFfRG::zaqbq7_qcd::tr67(fenv));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den11 + _den12 + _den13 + _den14 + _den15 + _den2 + _den3 + _den4 + _den5 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p2 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
      // clang-format off
      _acc += 0.005208333333333332 * fma(4., _den14 * _den2 * _interp11 * _interp17 * _interp20 * _interp21 * powr<-5>(l1) * powr<-4>(p), fma(4., _den14 * _den2 * _den9 * _interp11 * _interp15 * _interp21 * _interp22 * _interp23 * powr<-5>(l1) * powr<-4>(p), fma(4., _den2 * _den9 * _interp11 * _interp24 * _interp25 * _interp26 * powr<-5>(l1) * powr<-4>(p), fma(4., _den14 * _den2 * _interp11 * _interp27 * _interp28 * _interp29 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den14 * _den2 * _den9 * _interp11 * _interp23 * _interp26 * _interp29 * _interp30 * powr<-5>(l1) * powr<-4>(p), fma(2., _den2 * _den9 * _interp11 * _interp14 * _interp31 * _interp32 * powr<-5>(l1) * powr<-4>(p), fma(2., _den14 * _den2 * _den9 * _interp11 * _interp21 * _interp23 * _interp32 * _interp33 * powr<-5>(l1) * powr<-4>(p), fma(2., _den14 * _den2 * _interp11 * _interp20 * _interp34 * _interp35 * powr<-5>(l1) * powr<-4>(p), fma(2., _den14 * _den2 * _den9 * _interp11 * _interp15 * _interp23 * _interp35 * _interp36 * powr<-5>(l1) * powr<-4>(p), fma(_den14, _den2 * _den9 * _interp11 * _interp23 * _interp32 * _interp35 * _interp37 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp21 * _interp38 * _interp39 * _interp40 * powr<-5>(l1) * powr<-4>(p), fma(4., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp39 * _interp43 * _interp46 * _interp47 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp29 * _interp40 * _interp48 * _interp49 * powr<-5>(l1) * powr<-4>(p), fma(4., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp47 * _interp49 * _interp50 * _interp51 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den14 * _den2 * _interp10 * _interp17 * _interp20 * _interp21 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den14 * _den2 * _den9 * _interp10 * _interp15 * _interp21 * _interp22 * _interp23 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den2 * _den9 * _interp10 * _interp24 * _interp25 * _interp26 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den14 * _den2 * _interp10 * _interp27 * _interp28 * _interp29 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den14 * _den2 * _den9 * _interp10 * _interp23 * _interp26 * _interp29 * _interp30 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den2 * _den9 * _interp10 * _interp14 * _interp31 * _interp32 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den14 * _den2 * _den9 * _interp10 * _interp21 * _interp23 * _interp32 * _interp33 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den14 * _den2 * _interp10 * _interp20 * _interp34 * _interp35 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den14 * _den2 * _den9 * _interp10 * _interp15 * _interp23 * _interp35 * _interp36 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-1., _den14 * _den2 * _den9 * _interp10 * _interp23 * _interp32 * _interp35 * _interp37 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp21 * _interp38 * _interp39 * _interp40 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp39 * _interp43 * _interp46 * _interp47 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp29 * _interp40 * _interp48 * _interp49 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp47 * _interp49 * _interp50 * _interp51 * etaQ * powr<-5>(l1) * powr<-4>(p), 0.))))))))))))))))))))))))))));
      // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += 0.005208333333333332 * fma(-4., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp21 * _interp49 * _interp52 * _interp53 * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp29 * _interp39 * _interp53 * _interp54 * powr<-5>(l1) * powr<-4>(p), fma(4., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp46 * _interp49 * _interp55 * _interp56 * powr<-5>(l1) * powr<-4>(p), fma(-4., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp39 * _interp51 * _interp56 * _interp57 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp21 * _interp40 * _interp58 * _interp59 * powr<-5>(l1) * powr<-4>(p), fma(2., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp46 * _interp47 * _interp59 * _interp60 * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp29 * _interp53 * _interp59 * _interp61 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp51 * _interp56 * _interp59 * _interp62 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp35 * _interp39 * _interp40 * _interp63 * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp35 * _interp49 * _interp53 * _interp64 * powr<-5>(l1) * powr<-4>(p), fma(-1., _den11 * _den14 * _den2 * _den4 * _interp11 * _interp35 * _interp40 * _interp59 * _interp65 * powr<-5>(l1) * powr<-4>(p), fma(2., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp39 * _interp47 * _interp66 * _interp67 * powr<-5>(l1) * powr<-4>(p), fma(2., _den10 * _den11 * _den2 * _den4 * _interp11 * _interp49 * _interp56 * _interp67 * _interp68 * powr<-5>(l1) * powr<-4>(p), fma(_den10, _den11 * _den2 * _den4 * _interp11 * _interp47 * _interp59 * _interp67 * _interp69 * powr<-5>(l1) * powr<-4>(p), fma(4., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp21 * _interp49 * _interp52 * _interp53 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp29 * _interp39 * _interp53 * _interp54 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-4., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp46 * _interp49 * _interp55 * _interp56 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(4., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp39 * _interp51 * _interp56 * _interp57 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp21 * _interp40 * _interp58 * _interp59 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp46 * _interp47 * _interp59 * _interp60 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp29 * _interp53 * _interp59 * _interp61 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp51 * _interp56 * _interp59 * _interp62 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp35 * _interp39 * _interp40 * _interp63 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(2., _den11 * _den14 * _den2 * _den4 * _interp10 * _interp35 * _interp49 * _interp53 * _interp64 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(_den11, _den14 * _den2 * _den4 * _interp10 * _interp35 * _interp40 * _interp59 * _interp65 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp39 * _interp47 * _interp66 * _interp67 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-2., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp49 * _interp56 * _interp67 * _interp68 * etaQ * powr<-5>(l1) * powr<-4>(p), fma(-1., _den10 * _den11 * _den2 * _den4 * _interp10 * _interp47 * _interp59 * _interp67 * _interp69 * etaQ * powr<-5>(l1) * powr<-4>(p), 0.))))))))))))))))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(-0.01041666666666667, _den11 * _den14 * _den2 * _den4 * _interp21 * _interp39 * _interp70 * _interp71 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(-0.01041666666666667, _den11 * _den14 * _den2 * _den4 * _interp29 * _interp49 * _interp71 * _interp72 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(-0.005208333333333332, _den11 * _den14 * _den2 * _den4 * _interp21 * _interp59 * _interp71 * _interp73 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(-0.005208333333333332, _den11 * _den14 * _den2 * _den4 * _interp35 * _interp39 * _interp71 * _interp74 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(-0.002604166666666666, _den11 * _den14 * _den2 * _den4 * _interp35 * _interp59 * _interp71 * _interp75 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.01041666666666667, _den10 * _den11 * _den2 * _den4 * _interp39 * _interp46 * _interp76 * _interp77 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.01041666666666667, _den10 * _den11 * _den2 * _den4 * _interp49 * _interp51 * _interp77 * _interp78 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.005208333333333332, _den10 * _den11 * _den2 * _den4 * _interp46 * _interp59 * _interp77 * _interp79 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.005208333333333332, _den10 * _den11 * _den2 * _den4 * _interp39 * _interp67 * _interp77 * _interp80 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.002604166666666666, _den10 * _den11 * _den2 * _den4 * _interp59 * _interp67 * _interp77 * _interp81 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(0.02083333333333333, _den2 * _den9 * _interp14 * _interp15 * _interp9 * (_interp11 - _interp10 * etaQ) * powr<-5>(l1) * powr<-4>(p), fma(-0.01041666666666667, _den12 * _den15 * _den3 * _den5 * _den7 * _interp108 * _interp109 * _interp15 * _interp21 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den12 * _den15 * _den3 * _den5 * _den7 * _interp100 * _interp101 * _interp21 * _interp26 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den12 * _den15 * _den3 * _den5 * _den7 * _interp101 * _interp102 * _interp15 * _interp29 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 4
        const auto _den6 = powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        // clang-format off
        _acc += fma(-0.01041666666666667, _den12 * _den15 * _den3 * _den5 * _den7 * _interp109 * _interp110 * _interp26 * _interp29 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.005208333333333332, _den12 * _den15 * _den3 * _den5 * _den7 * _interp109 * _interp111 * _interp21 * _interp32 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den12 * _den15 * _den3 * _den5 * _den7 * _interp101 * _interp104 * _interp29 * _interp32 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.005208333333333332, _den12 * _den15 * _den3 * _den5 * _den7 * _interp109 * _interp112 * _interp15 * _interp35 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den12 * _den15 * _den3 * _den5 * _den7 * _interp101 * _interp106 * _interp26 * _interp35 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.002604166666666666, _den12 * _den15 * _den3 * _den5 * _den7 * _interp109 * _interp113 * _interp32 * _interp35 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den13 * _den3 * _den6 * _den8 * _interp114 * _interp46 * _interp47 * _interp7 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den13 * _den3 * _den6 * _den8 * _interp117 * _interp51 * _interp56 * _interp7 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den13 * _den3 * _den6 * _den8 * _interp118 * _interp47 * _interp67 * _interp7 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den13 * _den3 * _den6 * _den8 * _interp119 * _interp46 * _interp7 * _interp77 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.005208333333333332, _den13 * _den3 * _den6 * _den8 * _interp120 * _interp67 * _interp7 * _interp77 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den3 * _den8 * _interp1 * _interp7 * _interp8 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den12 * _den3 * _den5 * _interp14 * _interp15 * _interp82 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den12 * _den3 * _den5 * _interp25 * _interp26 * _interp85 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 5
      // clang-format off
      _acc += fma(0.01041666666666667, _den12 * _den3 * _den5 * _interp14 * _interp32 * _interp86 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den15 * _den3 * _den7 * _interp20 * _interp21 * _interp87 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den15 * _den3 * _den7 * _den8 * _interp21 * _interp40 * _interp7 * _interp90 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den15 * _den3 * _den7 * _interp28 * _interp29 * _interp91 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.02083333333333333, _den15 * _den3 * _den7 * _den8 * _interp29 * _interp53 * _interp7 * _interp92 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den15 * _den3 * _den7 * _interp20 * _interp35 * _interp93 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den15 * _den3 * _den7 * _den8 * _interp35 * _interp40 * _interp7 * _interp94 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.01041666666666667, _den15 * _den3 * _den7 * _den8 * _interp21 * _interp7 * _interp71 * _interp95 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(0.005208333333333332, _den15 * _den3 * _den7 * _den8 * _interp35 * _interp7 * _interp71 * _interp96 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den12 * _den15 * _den3 * _den5 * _den7 * _interp103 * _interp21 * _interp32 * _interp98 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.01041666666666667, _den12 * _den15 * _den3 * _den5 * _den7 * _interp105 * _interp15 * _interp35 * _interp98 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.005208333333333332, _den12 * _den15 * _den3 * _den5 * _den7 * _interp107 * _interp32 * _interp35 * _interp98 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den12 * _den15 * _den3 * _den5 * _den7 * _interp15 * _interp21 * _interp97 * _interp98 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), fma(-0.02083333333333333, _den12 * _den15 * _den3 * _den5 * _den7 * _interp26 * _interp29 * _interp98 * _interp99 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-4>(p), 0.))))))))))))));
      // clang-format on
      }
      return _acc;
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaqbq7_qcd::nenv) > 0 ? (DiFfRG::zaqbq7_qcd::nenv) : 1];
      DiFfRG::zaqbq7_qcd::fill(fenv, l1, cos1, cos2, p);
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix8 = ZA.index(l1);
      const auto _interp3 = ZA.at(_ix8);
      const auto _interp8 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp9 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _ix1 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1)))));
      const auto _interp12 = ZAAqbq2.at(_ix1);
      const auto _ix5 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp13 = ZAqbq1.at(_ix5);
      const auto _interp14 = zq.at(_ix8);
      const auto _interp15 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _ix10 = ZA.index(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp16 = ZA.at(_ix10);
      const auto _ix0 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp17 = ZAAqbq2.at(_ix0);
      const auto _ix2 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp18 = ZAqbq1.at(_ix2);
      const auto _interp20 = ZAAqbq1.at(_ix1);
      const auto _interp21 = ZAqbq4.at(_ix5);
      const auto _interp22 = ZAAqbq1.at(_ix0);
      const auto _interp23 = ZAqbq4.at(_ix2);
      const auto _interp24 = ZAqbq7.at(_ix5);
      const auto _interp25 = ZAqbq7.at(_ix2);
      const auto _ix4 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp26 = ZAqbq1.at(_ix4);
      const auto _ix6 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp27 = ZAqbq1.at(_ix6);
      const auto _interp28 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp29 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp30 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp31 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _ix7 = ZAqbq1.index(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = ZAqbq1.at(_ix7);
      const auto _ix3 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _interp33 = ZAqbq1.at(_ix3);
      const auto _interp34 = ZAqbq4.at(_ix4);
      const auto _interp35 = ZAqbq4.at(_ix7);
      const auto _interp36 = ZAqbq4.at(_ix6);
      const auto _interp37 = ZAqbq4.at(_ix3);
      const auto _interp38 = ZAqbq7.at(_ix4);
      const auto _interp39 = ZAqbq7.at(_ix7);
      const auto _interp40 = ZAqbq7.at(_ix6);
      const auto _interp41 = ZAqbq7.at(_ix3);
      const auto _interp42 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp43 = zq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp44 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp45 = zq.at(_ix10);
      const auto _interp48 = ZAqbq7(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.5, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp8, powr<-1>(l1), _interp14));
      const auto _den3 = powr<-2>(fma(_interp3, powr<2>(l1), fma(_interp1, nthk0, 0.)));
      const auto _den4 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den5 = powr<-1>(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den7 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den10 = powr<-1>(fma(_interp30, nthk0, fma(_interp31, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den11 = powr<-1>(fma(_interp28, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp29));
      const auto _den12 = powr<-1>(fma(_interp42, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))), _interp43));
      const auto _den14 = powr<-1>(fma(_interp15, nthk0, fma(_interp16, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den15 = powr<-1>(fma(_interp44, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))), _interp45));
      const auto _interp10 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp11 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _ix9 = ZA3.index(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp19 = ZA3.at(_ix9);
      const auto _interp46 = ZAqbq1.at(_ix9);
      const auto _interp47 = ZAqbq4.at(_ix9);
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp5 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp7 = ZAAqbq2(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
      const auto _interp49 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp50 = zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den13 = powr<-1>(_interp50 + _interp49 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(complex<double>(0.,-8.) + _den1 + _den10 + _den11 + _den12 + _den13 + _den14 + _den15 + _den2 + _den3 + _den4 + _den5 + _den7 + DiFfRG::zaqbq7_qcd::tr0(fenv) + DiFfRG::zaqbq7_qcd::tr10(fenv) + DiFfRG::zaqbq7_qcd::tr11(fenv) + DiFfRG::zaqbq7_qcd::tr12(fenv) + DiFfRG::zaqbq7_qcd::tr13(fenv) + DiFfRG::zaqbq7_qcd::tr14(fenv) + DiFfRG::zaqbq7_qcd::tr15(fenv) + DiFfRG::zaqbq7_qcd::tr16(fenv) + DiFfRG::zaqbq7_qcd::tr17(fenv) + DiFfRG::zaqbq7_qcd::tr18(fenv) + DiFfRG::zaqbq7_qcd::tr19(fenv) + DiFfRG::zaqbq7_qcd::tr1(fenv) + DiFfRG::zaqbq7_qcd::tr20(fenv) + DiFfRG::zaqbq7_qcd::tr21(fenv) + DiFfRG::zaqbq7_qcd::tr22(fenv) + DiFfRG::zaqbq7_qcd::tr23(fenv) + DiFfRG::zaqbq7_qcd::tr24(fenv) + DiFfRG::zaqbq7_qcd::tr25(fenv) + DiFfRG::zaqbq7_qcd::tr26(fenv) + DiFfRG::zaqbq7_qcd::tr27(fenv) + DiFfRG::zaqbq7_qcd::tr28(fenv) + DiFfRG::zaqbq7_qcd::tr29(fenv) + DiFfRG::zaqbq7_qcd::tr2(fenv) + DiFfRG::zaqbq7_qcd::tr30(fenv) + DiFfRG::zaqbq7_qcd::tr31(fenv) + DiFfRG::zaqbq7_qcd::tr32(fenv) + DiFfRG::zaqbq7_qcd::tr33(fenv) + DiFfRG::zaqbq7_qcd::tr34(fenv) + DiFfRG::zaqbq7_qcd::tr35(fenv) + DiFfRG::zaqbq7_qcd::tr36(fenv) + DiFfRG::zaqbq7_qcd::tr37(fenv) + DiFfRG::zaqbq7_qcd::tr38(fenv) + DiFfRG::zaqbq7_qcd::tr39(fenv) + DiFfRG::zaqbq7_qcd::tr3(fenv) + DiFfRG::zaqbq7_qcd::tr40(fenv) + DiFfRG::zaqbq7_qcd::tr41(fenv) + DiFfRG::zaqbq7_qcd::tr42(fenv) + DiFfRG::zaqbq7_qcd::tr43(fenv) + DiFfRG::zaqbq7_qcd::tr44(fenv) + DiFfRG::zaqbq7_qcd::tr45(fenv) + DiFfRG::zaqbq7_qcd::tr46(fenv) + DiFfRG::zaqbq7_qcd::tr47(fenv) + DiFfRG::zaqbq7_qcd::tr48(fenv) + DiFfRG::zaqbq7_qcd::tr49(fenv) + DiFfRG::zaqbq7_qcd::tr4(fenv) + DiFfRG::zaqbq7_qcd::tr50(fenv) + DiFfRG::zaqbq7_qcd::tr51(fenv) + DiFfRG::zaqbq7_qcd::tr52(fenv) + DiFfRG::zaqbq7_qcd::tr53(fenv) + DiFfRG::zaqbq7_qcd::tr54(fenv) + DiFfRG::zaqbq7_qcd::tr55(fenv) + DiFfRG::zaqbq7_qcd::tr56(fenv) + DiFfRG::zaqbq7_qcd::tr57(fenv) + DiFfRG::zaqbq7_qcd::tr58(fenv) + DiFfRG::zaqbq7_qcd::tr59(fenv) + DiFfRG::zaqbq7_qcd::tr5(fenv) + DiFfRG::zaqbq7_qcd::tr60(fenv) + DiFfRG::zaqbq7_qcd::tr61(fenv) + DiFfRG::zaqbq7_qcd::tr62(fenv) + DiFfRG::zaqbq7_qcd::tr63(fenv) + DiFfRG::zaqbq7_qcd::tr64(fenv) + DiFfRG::zaqbq7_qcd::tr65(fenv) + DiFfRG::zaqbq7_qcd::tr66(fenv) + DiFfRG::zaqbq7_qcd::tr67(fenv) + DiFfRG::zaqbq7_qcd::tr68(fenv) + DiFfRG::zaqbq7_qcd::tr69(fenv) + DiFfRG::zaqbq7_qcd::tr6(fenv) + DiFfRG::zaqbq7_qcd::tr7(fenv) + DiFfRG::zaqbq7_qcd::tr8(fenv) + DiFfRG::zaqbq7_qcd::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den9 = powr<-1>(_interp10 * nthk0 + _interp11 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _cse1_k1 = powr<-5>(l1);
        const auto _cse2_k1 = powr<-4>(p);
        const auto _cse3_k1 = -_interp8 * etaQ;
        const auto _cse4_k1 = _cse3_k1 + _interp9;
        // clang-format off
        _acc += 0.005208333333333332 * fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr1(fenv) * _interp12 * _interp13, fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * DiFfRG::zaqbq7_qcd::tr2(fenv) * _interp17 * _interp18, fma(complex<double>(0.,-4.), _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr3(fenv) * _interp13 * _interp18 * _interp19, fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr4(fenv) * _interp20 * _interp21, fma(complex<double>(0.,4.), _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr6(fenv) * _interp19 * _interp21 * _interp23, fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * DiFfRG::zaqbq7_qcd::tr5(fenv) * _interp22 * _interp23, fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr7(fenv) * _interp12 * _interp24, fma(complex<double>(0.,-2.), _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr8(fenv) * _interp18 * _interp19 * _interp24, fma(2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * DiFfRG::zaqbq7_qcd::tr9(fenv) * _interp17 * _interp25, fma(complex<double>(0.,-2.), _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr10(fenv) * _interp13 * _interp19 * _interp25, fma(complex<double>(0.,-1.), _cse1_k1 * _cse2_k1 * _cse4_k1 * _den14 * _den2 * _den9 * DiFfRG::zaqbq7_qcd::tr11(fenv) * _interp19 * _interp24 * _interp25, fma(-4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr13(fenv) * _interp18 * _interp26 * _interp27, fma(-2., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr35(fenv) * _interp25 * _interp26 * _interp27, fma(4., _cse1_k1 * _cse2_k1 * _cse4_k1 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr12(fenv) * _interp26 * _interp32 * _interp33, 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _cse1_k2 = powr<-5>(l1);
        const auto _cse2_k2 = powr<-4>(p);
        const auto _cse3_k2 = -_interp8 * etaQ;
        const auto _cse4_k2 = _cse3_k2 + _interp9;
        // clang-format off
        _acc += 0.005208333333333332 * fma(-4., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr18(fenv) * _interp23 * _interp27 * _interp34, fma(4., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr14(fenv) * _interp33 * _interp34 * _interp35, fma(4., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr19(fenv) * _interp23 * _interp26 * _interp36, fma(-4., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr15(fenv) * _interp18 * _interp34 * _interp36, fma(-2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr36(fenv) * _interp25 * _interp34 * _interp36, fma(4., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr16(fenv) * _interp32 * _interp34 * _interp37, fma(-4., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr17(fenv) * _interp26 * _interp35 * _interp37, fma(-2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr21(fenv) * _interp18 * _interp27 * _interp38, fma(-1., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr37(fenv) * _interp25 * _interp27 * _interp38, fma(2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr20(fenv) * _interp32 * _interp33 * _interp38, fma(2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr23(fenv) * _interp23 * _interp36 * _interp38, fma(-2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr22(fenv) * _interp35 * _interp37 * _interp38, fma(2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr24(fenv) * _interp26 * _interp33 * _interp39, fma(2., _cse1_k2 * _cse2_k2 * _cse4_k2 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr25(fenv) * _interp34 * _interp37 * _interp39, 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 3
        const auto _cse1_k3 = powr<-4>(p);
        const auto _cse2_k3 = _interp2 * nthk0;
        const auto _cse3_k3 = -nthk0;
        const auto _cse4_k3 = _cse3_k3 + nthk2;
        const auto _cse5_k3 = 50. * _cse4_k3 * _den1 * powr<6>(k);
        const auto _cse6_k3 = _cse5_k3 + nthk1;
        const auto _cse7_k3 = _cse6_k3 * _interp1;
        const auto _cse8_k3 = _cse2_k3 + _cse7_k3;
        const auto _cse9_k3 = powr<-5>(l1);
        const auto _cse10_k3 = -_interp8 * etaQ;
        const auto _cse11_k3 = _cse10_k3 + _interp9;
        // clang-format off
        _acc += 0.002604166666666666 * fma(8., _cse1_k3 * _cse8_k3 * _den12 * _den3 * _den5 * DiFfRG::zaqbq7_qcd::tr40(fenv) * _interp12 * _interp13, fma(8., _cse1_k3 * _cse8_k3 * _den15 * _den3 * _den7 * DiFfRG::zaqbq7_qcd::tr48(fenv) * _interp17 * _interp18, fma(8., _cse1_k3 * _cse8_k3 * _den12 * _den3 * _den5 * DiFfRG::zaqbq7_qcd::tr41(fenv) * _interp20 * _interp21, fma(2., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr26(fenv) * _interp33 * _interp38 * _interp39, fma(-4., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr27(fenv) * _interp18 * _interp26 * _interp40, fma(-2., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr38(fenv) * _interp25 * _interp26 * _interp40, fma(-4., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr28(fenv) * _interp23 * _interp34 * _interp40, fma(-2., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr29(fenv) * _interp18 * _interp38 * _interp40, fma(-1., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den11 * _den14 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr39(fenv) * _interp25 * _interp38 * _interp40, fma(4., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr30(fenv) * _interp26 * _interp32 * _interp41, fma(4., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr31(fenv) * _interp34 * _interp35 * _interp41, fma(2., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr32(fenv) * _interp32 * _interp38 * _interp41, fma(2., _cse11_k3 * _cse1_k3 * _cse9_k3 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr33(fenv) * _interp26 * _interp39 * _interp41, fma(_cse1_k3, _cse11_k3 * _cse9_k3 * _den10 * _den11 * _den2 * _den4 * DiFfRG::zaqbq7_qcd::tr34(fenv) * _interp38 * _interp39 * _interp41, 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 4
        const auto _cse1_k4 = powr<-4>(p);
        const auto _cse2_k4 = _interp2 * nthk0;
        const auto _cse3_k4 = -nthk0;
        const auto _cse4_k4 = _cse3_k4 + nthk2;
        const auto _cse5_k4 = 50. * _cse4_k4 * _den1 * powr<6>(k);
        const auto _cse6_k4 = _cse5_k4 + nthk1;
        const auto _cse7_k4 = _cse6_k4 * _interp1;
        const auto _cse8_k4 = _cse2_k4 + _cse7_k4;
        // clang-format off
        _acc += 0.005208333333333332 * fma(-4., _cse1_k4 * _cse8_k4 * _den15 * _den3 * _den7 * DiFfRG::zaqbq7_qcd::tr50(fenv) * _interp22 * _interp23, fma(2., _cse1_k4 * _cse8_k4 * _den12 * _den3 * _den5 * DiFfRG::zaqbq7_qcd::tr42(fenv) * _interp12 * _interp24, fma(2., _cse1_k4 * _cse8_k4 * _den15 * _den3 * _den7 * DiFfRG::zaqbq7_qcd::tr53(fenv) * _interp17 * _interp25, fma(-4., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr56(fenv) * _interp13 * _interp18 * _interp46, fma(-4., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr58(fenv) * _interp21 * _interp23 * _interp46, fma(-2., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr60(fenv) * _interp18 * _interp24 * _interp46, fma(-2., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr65(fenv) * _interp13 * _interp25 * _interp46, fma(-1., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr67(fenv) * _interp24 * _interp25 * _interp46, fma(-4., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr57(fenv) * _interp18 * _interp21 * _interp47, fma(4., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr59(fenv) * _interp13 * _interp23 * _interp47, fma(2., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr61(fenv) * _interp23 * _interp24 * _interp47, fma(-2., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr66(fenv) * _interp21 * _interp25 * _interp47, fma(-2., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr62(fenv) * _interp13 * _interp18 * _interp48, fma(-2., _cse1_k4 * _cse8_k4 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr63(fenv) * _interp21 * _interp23 * _interp48, 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 5
        const auto _den6 = powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den8 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _cse1_k5 = powr<-4>(p);
        const auto _cse2_k5 = _interp2 * nthk0;
        const auto _cse3_k5 = -nthk0;
        const auto _cse4_k5 = _cse3_k5 + nthk2;
        const auto _cse5_k5 = 50. * _cse4_k5 * _den1 * powr<6>(k);
        // clang-format off
        _acc += 0.002604166666666666 * fma(-2., _cse1_k5 * _cse2_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr64(fenv) * _interp18 * _interp24 * _interp48, fma(-2., _cse1_k5 * _cse5_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr64(fenv) * _interp1 * _interp18 * _interp24 * _interp48, fma(-2., _cse1_k5 * _cse2_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr68(fenv) * _interp13 * _interp25 * _interp48, fma(-2., _cse1_k5 * _cse5_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr68(fenv) * _interp1 * _interp13 * _interp25 * _interp48, fma(-1., _cse1_k5 * _cse2_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr69(fenv) * _interp24 * _interp25 * _interp48, fma(-1., _cse1_k5 * _cse5_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr69(fenv) * _interp1 * _interp24 * _interp25 * _interp48, fma(complex<double>(0.,-8.), _cse1_k5 * _cse2_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr49(fenv) * _interp18 * _interp27 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse5_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr49(fenv) * _interp1 * _interp18 * _interp27 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse2_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr54(fenv) * _interp25 * _interp27 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse5_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr54(fenv) * _interp1 * _interp25 * _interp27 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse2_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr43(fenv) * _interp32 * _interp33 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse5_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr43(fenv) * _interp1 * _interp32 * _interp33 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse2_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr51(fenv) * _interp23 * _interp36 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse5_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr51(fenv) * _interp1 * _interp23 * _interp36 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse2_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr44(fenv) * _interp35 * _interp37 * _interp6, fma(complex<double>(0.,-8.), _cse1_k5 * _cse5_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr44(fenv) * _interp1 * _interp35 * _interp37 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse2_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr45(fenv) * _interp33 * _interp39 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse5_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr45(fenv) * _interp1 * _interp33 * _interp39 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse2_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr52(fenv) * _interp18 * _interp40 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse5_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr52(fenv) * _interp1 * _interp18 * _interp40 * _interp6, fma(complex<double>(0.,-2.), _cse1_k5 * _cse2_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr55(fenv) * _interp25 * _interp40 * _interp6, fma(complex<double>(0.,-2.), _cse1_k5 * _cse5_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr55(fenv) * _interp1 * _interp25 * _interp40 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse2_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr46(fenv) * _interp32 * _interp41 * _interp6, fma(complex<double>(0.,-4.), _cse1_k5 * _cse5_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr46(fenv) * _interp1 * _interp32 * _interp41 * _interp6, fma(complex<double>(0.,-2.), _cse1_k5 * _cse2_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr47(fenv) * _interp39 * _interp41 * _interp6, fma(complex<double>(0.,-2.), _cse1_k5 * _cse5_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr47(fenv) * _interp1 * _interp39 * _interp41 * _interp6, fma(complex<double>(0.,8.), _cse1_k5 * _cse2_k5 * _den3 * _den8 * DiFfRG::zaqbq7_qcd::tr0(fenv) * _interp6 * _interp7, fma(complex<double>(0.,8.), _cse1_k5 * _cse5_k5 * _den3 * _den8 * DiFfRG::zaqbq7_qcd::tr0(fenv) * _interp1 * _interp6 * _interp7, fma(-2., _cse1_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr64(fenv) * _interp1 * _interp18 * _interp24 * _interp48 * nthk1, fma(-2., _cse1_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr68(fenv) * _interp1 * _interp13 * _interp25 * _interp48 * nthk1, fma(-1., _cse1_k5 * _den12 * _den15 * _den3 * _den5 * _den7 * DiFfRG::zaqbq7_qcd::tr69(fenv) * _interp1 * _interp24 * _interp25 * _interp48 * nthk1, fma(complex<double>(0.,-8.), _cse1_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr49(fenv) * _interp1 * _interp18 * _interp27 * _interp6 * nthk1, fma(complex<double>(0.,-4.), _cse1_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr54(fenv) * _interp1 * _interp25 * _interp27 * _interp6 * nthk1, fma(complex<double>(0.,-8.), _cse1_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr43(fenv) * _interp1 * _interp32 * _interp33 * _interp6 * nthk1, fma(complex<double>(0.,-8.), _cse1_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr51(fenv) * _interp1 * _interp23 * _interp36 * _interp6 * nthk1, fma(complex<double>(0.,-8.), _cse1_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr44(fenv) * _interp1 * _interp35 * _interp37 * _interp6 * nthk1, fma(complex<double>(0.,-4.), _cse1_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr45(fenv) * _interp1 * _interp33 * _interp39 * _interp6 * nthk1, fma(complex<double>(0.,-4.), _cse1_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr52(fenv) * _interp1 * _interp18 * _interp40 * _interp6 * nthk1, fma(complex<double>(0.,-2.), _cse1_k5 * _den15 * _den3 * _den7 * _den8 * DiFfRG::zaqbq7_qcd::tr55(fenv) * _interp1 * _interp25 * _interp40 * _interp6 * nthk1, fma(complex<double>(0.,-4.), _cse1_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr46(fenv) * _interp1 * _interp32 * _interp41 * _interp6 * nthk1, fma(complex<double>(0.,-2.), _cse1_k5 * _den13 * _den3 * _den6 * _den8 * DiFfRG::zaqbq7_qcd::tr47(fenv) * _interp1 * _interp39 * _interp41 * _interp6 * nthk1, fma(complex<double>(0.,8.), _cse1_k5 * _den3 * _den8 * DiFfRG::zaqbq7_qcd::tr0(fenv) * _interp1 * _interp6 * _interp7 * nthk1, 0.))))))))))))))))))))))))))))))))))))))))));
        // clang-format on
      }
      return _acc;
    }
    #endif

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const auto _interp1 = ZAqbq7(p);
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

    static KOKKOS_INLINE_FUNCTION double ntRe(double x) { return x; }
    template <class T> static KOKKOS_INLINE_FUNCTION double ntRe(const T &z) { return z.real(); }
    static KOKKOS_INLINE_FUNCTION double ntIm(double) { return 0.0; }
    template <class T> static KOKKOS_INLINE_FUNCTION double ntIm(const T &z) { return z.imag(); }
  };
}
using DiFfRG::ZAqbq7_kernel;