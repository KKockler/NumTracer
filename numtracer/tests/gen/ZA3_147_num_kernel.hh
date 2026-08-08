#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "ZA3_147_num_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZA3_147_num_kernel
  {
    public:
    #if NT_ZA3_147_NUM_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
      const double dr_0 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_1 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_2 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_3 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_4 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_5 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_6 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_7 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_8 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = Mq(l1);
      const double dr_14 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      const double dr_15 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_16 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_17 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_20 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_21 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_22 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_23 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_24 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_25 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_26 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_27 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_28 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_29 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_30 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23, dr_24, dr_25, dr_26, dr_27, dr_28, dr_29, dr_30);
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp9 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp10 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp11 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp16 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp17 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp18 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp3, _interp4, fma(_interp7, powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp16, _interp3, fma(_interp17, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp3, _interp8, fma(_interp9, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den11 = powr<-1>(fma(_interp10, _interp3, fma(_interp11, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _interp1 = ntRe(DiFfRG::za3_147_num::tr1(fenv));
      const auto _interp14 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp15 = ntRe(DiFfRG::za3_147_num::tr2(fenv));
      const auto _interp19 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = ntRe(DiFfRG::za3_147_num::tr0(fenv));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp23 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp24 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp31 = ntRe(DiFfRG::za3_147_num::tr8(fenv));
      const auto _interp32 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp33 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp34 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp35 = Zc(k);
      const auto _interp36 = dtZc(k);
      const auto _interp37 = Zc(1.02 * k);
      const auto _interp38 = Zc(l1);
      const auto _interp39 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp40 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp41 = ntRe(DiFfRG::za3_147_num::tr9(fenv));
      const auto _interp42 = ntRe(DiFfRG::za3_147_num::tr10(fenv));
      const auto _interp43 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp44 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp45 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp46 = ntRe(DiFfRG::za3_147_num::tr11(fenv));
      const auto _interp47 = ntRe(DiFfRG::za3_147_num::tr6(fenv));
      const auto _interp48 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp49 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp50 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp51 = ntRe(DiFfRG::za3_147_num::tr7(fenv));
      const auto _interp25 = ntRe(DiFfRG::za3_147_num::tr3(fenv));
      const auto _interp26 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
      const auto _interp27 = ntRe(DiFfRG::za3_147_num::tr4(fenv));
      const auto _interp28 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p2 * l1 * p + 3. * powr<2>(p)));
      const auto _interp29 = ntRe(DiFfRG::za3_147_num::tr5(fenv));
      const auto _interp30 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp52 = ntRe(DiFfRG::za3_147_num::tr14(fenv));
      const auto _interp53 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp54 = Zq(k);
      const auto _interp55 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp56 = dtZq(k);
      const auto _interp57 = Zq(1.02 * k);
      const auto _interp58 = Mq(l1);
      const auto _interp59 = Zq(l1);
      const auto _interp60 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp61 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp62 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp63 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp64 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp65 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp66 = ntRe(DiFfRG::za3_147_num::tr15(fenv));
      const auto _interp67 = ntRe(DiFfRG::za3_147_num::tr16(fenv));
      const auto _interp68 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp69 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp70 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp71 = ntRe(DiFfRG::za3_147_num::tr17(fenv));
      const auto _interp72 = ntRe(DiFfRG::za3_147_num::tr12(fenv));
      const auto _interp73 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp74 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp75 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp76 = ntRe(DiFfRG::za3_147_num::tr13(fenv));
      const auto _den13 = powr<-1>(powr<2>(_interp60) + powr<2>(_interp54 * _interp61 + _interp62 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den14 = powr<-1>(powr<2>(_interp68) + powr<2>(_interp54 * _interp69 + _interp70 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den15 = powr<-1>(powr<2>(_interp73) + powr<2>(_interp54 * _interp74 + _interp75 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den16 = powr<-1>(powr<2>(_interp63) + powr<2>(_interp65) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp54 * _interp64 * (_interp54 * _interp64 + 2. * _interp65 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den11 + _den13 + _den14 + _den15 + _den16 + _den2 + _den4 + _den5 + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den3 = powr<-2>(_interp35 * _interp4 + _interp38 * powr<2>(l1));
        const auto _den6 = powr<-1>(_interp21 * _interp3 + _interp22 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den7 = -powr<-1>(_interp16 * _interp35 + _interp39 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den8 = -powr<-1>(_interp35 * _interp8 + _interp45 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den9 = -powr<-1>(_interp21 * _interp35 + _interp50 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp10 * _interp35 + _interp40 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1 = powr<-2>(p);
        const auto _cse2 = -_interp3;
        const auto _cse3 = _cse2 + _interp6;
        const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
        const auto _cse5 = _cse4 + _interp5;
        const auto _cse6 = _cse5 * _interp4;
        const auto _cse7 = _interp2 * _interp3;
        const auto _cse8 = _cse6 + _cse7;
        // clang-format off
        _acc += fma(-0.001683501683501684, _cse1 * _cse8 * _den11 * _den2 * _den5 * _interp1 * _interp12 * _interp13 * _interp14, fma(-0.001683501683501684, _cse1 * _cse8 * _den11 * _den2 * _den4 * _interp13 * _interp15 * _interp18 * _interp19, fma(-0.001683501683501684, _cse1 * _cse8 * _den2 * _den4 * _den6 * _interp18 * _interp20 * _interp23 * _interp24, fma(-0.001683501683501684, _cse1 * _den12 * _den3 * _den7 * _interp31 * _interp32 * _interp33 * _interp34 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4), fma(0.001683501683501684, _cse1 * _den12 * _den3 * _den7 * _interp32 * _interp33 * _interp34 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp41, fma(-0.001683501683501684, _cse1 * _den12 * _den3 * _den8 * _interp33 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp42 * _interp43 * _interp44, fma(0.001683501683501684, _cse1 * _den12 * _den3 * _den8 * _interp33 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp43 * _interp44 * _interp46, fma(0.001683501683501684, _cse1 * _den3 * _den7 * _den9 * _interp32 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp47 * _interp48 * _interp49, fma(-0.001683501683501684, _cse1 * _den3 * _den7 * _den9 * _interp32 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp48 * _interp49 * _interp51, 0.)))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den10 = powr<-2>(powr<2>(_interp58) + powr<2>(_interp54 * _interp55 + _interp59 * l1));
        const auto _cse1 = powr<-1>(l1);
        const auto _cse2 = powr<-2>(p);
        const auto _cse3 = -_cse1 * _interp53 * _interp54;
        const auto _cse4 = -_interp54;
        const auto _cse5 = _cse4 + _interp57;
        const auto _cse6 = 50. * _cse5;
        const auto _cse7 = _cse6 + _interp56;
        const auto _cse8 = -_cse1 * _cse7 * _interp55;
        // clang-format off
        _acc += fma(0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den16 * _interp52, fma(-0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den16 * _interp66, fma(0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den14 * _den16 * _interp67, fma(-0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den14 * _den16 * _interp71, fma(-0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den15 * _interp72, fma(0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den15 * _interp76, fma(-0.001683501683501684, _cse2 * _den2 * _den4 * _interp18 * _interp25 * _interp26 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(-0.001683501683501684, _cse2 * _den2 * _den5 * _interp12 * _interp27 * _interp28 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(-0.001683501683501684, _cse2 * _den11 * _den2 * _interp13 * _interp29 * _interp30 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.)))))))));
        // clang-format on
      }
      return _acc;
    }
    #elif NT_ZA3_147_NUM_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
      const double dr_0 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_1 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_2 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_3 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_4 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_5 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_6 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_7 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_8 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = Mq(l1);
      const double dr_14 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      const double dr_15 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_16 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_17 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_20 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_21 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_22 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_23 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_24 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_25 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_26 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_27 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_28 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_29 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_30 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23, dr_24, dr_25, dr_26, dr_27, dr_28, dr_29, dr_30);
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp9 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp10 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp11 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp13 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp16 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp17 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp18 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp3, _interp4, fma(_interp7, powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp16, _interp3, fma(_interp17, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp3, _interp8, fma(_interp9, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den11 = powr<-1>(fma(_interp10, _interp3, fma(_interp11, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _interp1 = ntRe(DiFfRG::za3_147_num::tr1(fenv));
      const auto _interp14 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp15 = ntRe(DiFfRG::za3_147_num::tr2(fenv));
      const auto _interp19 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp20 = ntRe(DiFfRG::za3_147_num::tr0(fenv));
      const auto _interp21 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp22 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp23 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp24 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp31 = ntRe(DiFfRG::za3_147_num::tr8(fenv));
      const auto _interp32 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp33 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp34 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp35 = Zc(k);
      const auto _interp36 = dtZc(k);
      const auto _interp37 = Zc(1.02 * k);
      const auto _interp38 = Zc(l1);
      const auto _interp39 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp40 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp41 = ntRe(DiFfRG::za3_147_num::tr9(fenv));
      const auto _interp42 = ntRe(DiFfRG::za3_147_num::tr10(fenv));
      const auto _interp43 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp44 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp45 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp46 = ntRe(DiFfRG::za3_147_num::tr11(fenv));
      const auto _interp47 = ntRe(DiFfRG::za3_147_num::tr6(fenv));
      const auto _interp48 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp49 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp50 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp51 = ntRe(DiFfRG::za3_147_num::tr7(fenv));
      const auto _interp25 = ntRe(DiFfRG::za3_147_num::tr3(fenv));
      const auto _interp26 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
      const auto _interp27 = ntRe(DiFfRG::za3_147_num::tr4(fenv));
      const auto _interp28 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p2 * l1 * p + 3. * powr<2>(p)));
      const auto _interp29 = ntRe(DiFfRG::za3_147_num::tr5(fenv));
      const auto _interp30 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp52 = ntRe(DiFfRG::za3_147_num::tr14(fenv));
      const auto _interp53 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp54 = Zq(k);
      const auto _interp55 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp56 = dtZq(k);
      const auto _interp57 = Zq(1.02 * k);
      const auto _interp58 = Mq(l1);
      const auto _interp59 = Zq(l1);
      const auto _interp60 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp61 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp62 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp63 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp64 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp65 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp66 = ntRe(DiFfRG::za3_147_num::tr15(fenv));
      const auto _interp67 = ntRe(DiFfRG::za3_147_num::tr16(fenv));
      const auto _interp68 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp69 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp70 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp71 = ntRe(DiFfRG::za3_147_num::tr17(fenv));
      const auto _interp72 = ntRe(DiFfRG::za3_147_num::tr12(fenv));
      const auto _interp73 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp74 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp75 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp76 = ntRe(DiFfRG::za3_147_num::tr13(fenv));
      const auto _den13 = powr<-1>(powr<2>(_interp60) + powr<2>(_interp54 * _interp61 + _interp62 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den14 = powr<-1>(powr<2>(_interp68) + powr<2>(_interp54 * _interp69 + _interp70 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den15 = powr<-1>(powr<2>(_interp73) + powr<2>(_interp54 * _interp74 + _interp75 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den16 = powr<-1>(powr<2>(_interp63) + powr<2>(_interp65) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp54 * _interp64 * (_interp54 * _interp64 + 2. * _interp65 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den11 + _den13 + _den14 + _den15 + _den16 + _den2 + _den4 + _den5 + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den3 = powr<-2>(_interp35 * _interp4 + _interp38 * powr<2>(l1));
        const auto _den6 = powr<-1>(_interp21 * _interp3 + _interp22 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den7 = -powr<-1>(_interp16 * _interp35 + _interp39 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den8 = -powr<-1>(_interp35 * _interp8 + _interp45 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den9 = -powr<-1>(_interp21 * _interp35 + _interp50 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp10 * _interp35 + _interp40 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1 = powr<-2>(p);
        const auto _cse2 = -_interp3;
        const auto _cse3 = _cse2 + _interp6;
        const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
        const auto _cse5 = _cse4 + _interp5;
        const auto _cse6 = _cse5 * _interp4;
        const auto _cse7 = _interp2 * _interp3;
        const auto _cse8 = _cse6 + _cse7;
        // clang-format off
        _acc += fma(-0.001683501683501684, _cse1 * _cse8 * _den11 * _den2 * _den5 * _interp1 * _interp12 * _interp13 * _interp14, fma(-0.001683501683501684, _cse1 * _cse8 * _den11 * _den2 * _den4 * _interp13 * _interp15 * _interp18 * _interp19, fma(-0.001683501683501684, _cse1 * _cse8 * _den2 * _den4 * _den6 * _interp18 * _interp20 * _interp23 * _interp24, fma(-0.001683501683501684, _cse1 * _den12 * _den3 * _den7 * _interp31 * _interp32 * _interp33 * _interp34 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4), fma(0.001683501683501684, _cse1 * _den12 * _den3 * _den7 * _interp32 * _interp33 * _interp34 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp41, fma(-0.001683501683501684, _cse1 * _den12 * _den3 * _den8 * _interp33 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp42 * _interp43 * _interp44, fma(0.001683501683501684, _cse1 * _den12 * _den3 * _den8 * _interp33 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp43 * _interp44 * _interp46, fma(0.001683501683501684, _cse1 * _den3 * _den7 * _den9 * _interp32 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp47 * _interp48 * _interp49, fma(-0.001683501683501684, _cse1 * _den3 * _den7 * _den9 * _interp32 * (_interp2 * _interp35 + (_interp36 + 50. * (-_interp35 + _interp37)) * _interp4) * _interp48 * _interp49 * _interp51, 0.)))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den10 = powr<-2>(powr<2>(_interp58) + powr<2>(_interp54 * _interp55 + _interp59 * l1));
        const auto _cse1 = powr<-1>(l1);
        const auto _cse2 = powr<-2>(p);
        const auto _cse3 = -_cse1 * _interp53 * _interp54;
        const auto _cse4 = -_interp54;
        const auto _cse5 = _cse4 + _interp57;
        const auto _cse6 = 50. * _cse5;
        const auto _cse7 = _cse6 + _interp56;
        const auto _cse8 = -_cse1 * _cse7 * _interp55;
        // clang-format off
        _acc += fma(0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den16 * _interp52, fma(-0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den16 * _interp66, fma(0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den14 * _den16 * _interp67, fma(-0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den14 * _den16 * _interp71, fma(-0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den15 * _interp72, fma(0.003367003367003367, _cse2 * (_cse3 + _cse8) * _den10 * _den13 * _den15 * _interp76, fma(-0.001683501683501684, _cse2 * _den2 * _den4 * _interp18 * _interp25 * _interp26 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(-0.001683501683501684, _cse2 * _den2 * _den5 * _interp12 * _interp27 * _interp28 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(-0.001683501683501684, _cse2 * _den11 * _den2 * _interp13 * _interp29 * _interp30 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.)))))))));
        // clang-format on
      }
      return _acc;
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      double fenv[(DiFfRG::za3_147_num::nenv) > 0 ? (DiFfRG::za3_147_num::nenv) : 1];
      const double dr_0 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_1 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_2 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_3 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_4 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_5 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_6 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_7 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_8 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_10 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_13 = Mq(l1);
      const double dr_14 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      const double dr_15 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_16 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_17 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_20 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_21 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_22 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_23 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_24 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_25 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_26 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_27 = ZAqbq4(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_28 = ZAqbq7(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.5 * powr<2>(p)));
      const double dr_29 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_30 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      DiFfRG::za3_147_num::fill(fenv, l1, cos1, cos2, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23, dr_24, dr_25, dr_26, dr_27, dr_28, dr_29, dr_30);
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp7 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp8 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp9 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp10 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp11 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp12 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp14 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp15 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp16 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp2, _interp3, fma(_interp6, powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp14, _interp2, fma(_interp15, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp2, _interp7, fma(_interp8, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den11 = powr<-1>(fma(_interp2, _interp9, fma(_interp10, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _interp13 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp17 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp20 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp21 = ZA3(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp25 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp26 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp27 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp28 = Zc(k);
      const auto _interp29 = dtZc(k);
      const auto _interp30 = Zc(1.02 * k);
      const auto _interp31 = Zc(l1);
      const auto _interp32 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp33 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp34 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp35 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp36 = Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp37 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp38 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (-cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp39 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp22 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p1 * l1 * p + 3. * powr<2>(p)));
      const auto _interp23 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * cosl1p2 * l1 * p + 3. * powr<2>(p)));
      const auto _interp24 = ZA4(0.5 * sqrt(2. * powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 3. * powr<2>(p)));
      const auto _interp40 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp41 = Zq(k);
      const auto _interp42 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp43 = dtZq(k);
      const auto _interp44 = Zq(1.02 * k);
      const auto _interp45 = Mq(l1);
      const auto _interp46 = Zq(l1);
      const auto _interp47 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp48 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp49 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp50 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp51 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p));
      const auto _interp52 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
      const auto _interp53 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp54 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp55 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp56 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp57 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp58 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den13 = powr<-1>(powr<2>(_interp47) + powr<2>(_interp41 * _interp48 + _interp49 * sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den14 = powr<-1>(powr<2>(_interp53) + powr<2>(_interp41 * _interp54 + _interp55 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den15 = powr<-1>(powr<2>(_interp56) + powr<2>(_interp41 * _interp57 + _interp58 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den16 = powr<-1>(powr<2>(_interp50) + powr<2>(_interp52) * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + _interp41 * _interp51 * (_interp41 * _interp51 + 2. * _interp52 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den11 + _den13 + _den14 + _den15 + _den16 + _den2 + _den4 + _den5 + DiFfRG::za3_147_num::tr0(fenv) + DiFfRG::za3_147_num::tr10(fenv) + DiFfRG::za3_147_num::tr11(fenv) + DiFfRG::za3_147_num::tr12(fenv) + DiFfRG::za3_147_num::tr13(fenv) + DiFfRG::za3_147_num::tr14(fenv) + DiFfRG::za3_147_num::tr15(fenv) + DiFfRG::za3_147_num::tr16(fenv) + DiFfRG::za3_147_num::tr17(fenv) + DiFfRG::za3_147_num::tr1(fenv) + DiFfRG::za3_147_num::tr2(fenv) + DiFfRG::za3_147_num::tr3(fenv) + DiFfRG::za3_147_num::tr4(fenv) + DiFfRG::za3_147_num::tr5(fenv) + DiFfRG::za3_147_num::tr6(fenv) + DiFfRG::za3_147_num::tr7(fenv) + DiFfRG::za3_147_num::tr8(fenv) + DiFfRG::za3_147_num::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp6 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den3 = powr<-2>(_interp28 * _interp3 + _interp31 * powr<2>(l1));
        const auto _den6 = powr<-1>(_interp18 * _interp2 + _interp19 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den7 = -powr<-1>(_interp14 * _interp28 + _interp32 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den8 = -powr<-1>(_interp28 * _interp7 + _interp36 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den9 = -powr<-1>(_interp18 * _interp28 + _interp39 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den12 = -powr<-1>(_interp28 * _interp9 + _interp33 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1 = powr<-2>(p);
        const auto _cse2 = -_interp2;
        const auto _cse3 = _cse2 + _interp5;
        const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
        const auto _cse5 = _cse4 + _interp4;
        const auto _cse6 = _cse5 * _interp3;
        const auto _cse7 = _interp1 * _interp2;
        const auto _cse8 = _cse6 + _cse7;
        const auto _cse9 = -_interp28;
        const auto _cse10 = _cse9 + _interp30;
        const auto _cse11 = 50. * _cse10;
        const auto _cse12 = _cse11 + _interp29;
        const auto _cse13 = _cse12 * _interp3;
        const auto _cse14 = _interp1 * _interp28;
        const auto _cse15 = _cse13 + _cse14;
        // clang-format off
        _acc += 0.001683501683501684 * fma(-1., _cse1 * _cse8 * _den11 * _den2 * _den5 * DiFfRG::za3_147_num::tr1(fenv) * _interp11 * _interp12 * _interp13, fma(-1., _cse1 * _cse8 * _den11 * _den2 * _den4 * DiFfRG::za3_147_num::tr2(fenv) * _interp12 * _interp16 * _interp17, fma(-1., _cse1 * _cse8 * _den2 * _den4 * _den6 * DiFfRG::za3_147_num::tr0(fenv) * _interp16 * _interp20 * _interp21, fma(-1., _cse1 * _cse15 * _den12 * _den3 * _den7 * DiFfRG::za3_147_num::tr8(fenv) * _interp25 * _interp26 * _interp27, fma(_cse1, _cse15 * _den12 * _den3 * _den7 * DiFfRG::za3_147_num::tr9(fenv) * _interp25 * _interp26 * _interp27, fma(-1., _cse1 * _cse15 * _den12 * _den3 * _den8 * DiFfRG::za3_147_num::tr10(fenv) * _interp26 * _interp34 * _interp35, fma(_cse1, _cse15 * _den12 * _den3 * _den8 * DiFfRG::za3_147_num::tr11(fenv) * _interp26 * _interp34 * _interp35, fma(_cse1, _cse15 * _den3 * _den7 * _den9 * DiFfRG::za3_147_num::tr6(fenv) * _interp25 * _interp37 * _interp38, fma(-1., _cse1 * _cse15 * _den3 * _den7 * _den9 * DiFfRG::za3_147_num::tr7(fenv) * _interp25 * _interp37 * _interp38, 0.)))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den10 = powr<-2>(powr<2>(_interp45) + powr<2>(_interp41 * _interp42 + _interp46 * l1));
        const auto _cse1 = powr<-1>(l1);
        const auto _cse2 = powr<-2>(p);
        const auto _cse3 = -_cse1 * _interp40 * _interp41;
        const auto _cse4 = -_interp41;
        const auto _cse5 = _cse4 + _interp44;
        const auto _cse6 = 50. * _cse5;
        const auto _cse7 = _cse6 + _interp43;
        const auto _cse8 = -_cse1 * _cse7 * _interp42;
        const auto _cse9 = _cse3 + _cse8;
        const auto _cse10 = _interp1 * _interp2;
        const auto _cse11 = -_interp2;
        const auto _cse12 = _cse11 + _interp5;
        const auto _cse13 = 50. * _cse12 * _den1 * powr<6>(k);
        const auto _cse14 = _cse13 + _interp4;
        const auto _cse15 = _cse14 * _interp3;
        const auto _cse16 = _cse10 + _cse15;
        _acc += 0.001683501683501684 * fma(-2., _cse2 * _cse9 * _den10 * _den13 * _den15 * DiFfRG::za3_147_num::tr12(fenv), fma(2., _cse2 * _cse9 * _den10 * _den13 * _den15 * DiFfRG::za3_147_num::tr13(fenv), fma(2., _cse2 * _cse9 * _den10 * _den13 * _den16 * DiFfRG::za3_147_num::tr14(fenv), fma(-2., _cse2 * _cse9 * _den10 * _den13 * _den16 * DiFfRG::za3_147_num::tr15(fenv), fma(2., _cse2 * _cse9 * _den10 * _den14 * _den16 * DiFfRG::za3_147_num::tr16(fenv), fma(-2., _cse2 * _cse9 * _den10 * _den14 * _den16 * DiFfRG::za3_147_num::tr17(fenv), fma(-1., _cse16 * _cse2 * _den2 * _den4 * DiFfRG::za3_147_num::tr3(fenv) * _interp16 * _interp22, fma(-1., _cse16 * _cse2 * _den2 * _den5 * DiFfRG::za3_147_num::tr4(fenv) * _interp11 * _interp23, fma(-1., _cse16 * _cse2 * _den11 * _den2 * DiFfRG::za3_147_num::tr5(fenv) * _interp12 * _interp24, 0.)))))))));
      }
      return _acc;
    }
    #endif

    // clang-format off
    static inline auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      return 0.;
    }
    private:
    static inline double ntRe(double x) { return x; }
    template <class T> static inline double ntRe(const T &z) { return z.real(); }
    static inline double ntIm(double) { return 0.0; }
    template <class T> static inline double ntIm(const T &z) { return z.imag(); }
  };
}
using DiFfRG::ZA3_147_num_kernel;