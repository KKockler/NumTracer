#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZAAqbq2_kernel
  {
    public:
    using Regulator = REG;

    #if NT_ZAAQBQ2_QCD_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaaqbq2_qcd::nenv) > 0 ? (DiFfRG::zaaqbq2_qcd::nenv) : 1];
      DiFfRG::zaaqbq2_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix10 = ZA.index(l1);
      const auto _interp4 = ZA.at(_ix10);
      const auto _interp7 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp8 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _ix4 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp11 = ZAAqbq2.at(_ix4);
      const auto _interp13 = ZA4.at(_ix4);
      const auto _ix5 = ZAAqbq1.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp24 = ZAAqbq1.at(_ix5);
      const auto _interp27 = ZAAqbq2.at(_ix5);
      const auto _interp37 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp38 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp46 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp53 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp63 = ZAqbq1(0.816496580927726 * sqrt(fma(2., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp67 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp80 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den3 = powr<-2>(fma(_interp4, powr<2>(l1), fma(_interp2, nthk0, 0.)));
      const auto _den24 = powr<-1>(fma(_interp7, nthk0, fma(_interp8, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _interp12 = ntRe(DiFfRG::zaaqbq2_qcd::tr2(fenv));
      const auto _interp19 = ntRe(DiFfRG::zaaqbq2_qcd::tr3(fenv));
      const auto _interp20 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp21 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp22 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix14 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA.at(_ix14);
      const auto _interp25 = zq.at(_ix10);
      const auto _interp26 = ntRe(DiFfRG::zaaqbq2_qcd::tr4(fenv));
      const auto _interp28 = ntRe(DiFfRG::zaaqbq2_qcd::tr7(fenv));
      const auto _interp29 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _ix18 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp30 = ZA.at(_ix18);
      const auto _interp31 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _ix11 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp32 = ZA.at(_ix11);
      const auto _interp33 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix12 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp34 = ZA.at(_ix12);
      const auto _ix3 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp35 = ZA3.at(_ix3);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp36 = ZA3.at(_ix0);
      const auto _interp39 = ntRe(DiFfRG::zaaqbq2_qcd::tr8(fenv));
      const auto _ix6 = ZA4.index(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp40 = ZA4.at(_ix6);
      const auto _interp41 = ntRe(DiFfRG::zaaqbq2_qcd::tr9(fenv));
      const auto _interp42 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _ix19 = ZA.index(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp43 = ZA.at(_ix19);
      const auto _interp44 = ZAAqbq2(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _ix7 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp45 = ZAqbq1.at(_ix7);
      const auto _interp47 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix13 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp48 = zq.at(_ix13);
      const auto _interp49 = ntRe(DiFfRG::zaaqbq2_qcd::tr10(fenv));
      const auto _interp50 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix15 = ZA.index(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp51 = ZA.at(_ix15);
      const auto _interp52 = ZAAqbq2(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp54 = ntRe(DiFfRG::zaaqbq2_qcd::tr13(fenv));
      const auto _interp55 = ZAAqbq2(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _ix8 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp56 = ZAqbq1.at(_ix8);
      const auto _interp57 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix16 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp58 = zq.at(_ix16);
      const auto _interp59 = ntRe(DiFfRG::zaaqbq2_qcd::tr12(fenv));
      const auto _interp60 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp61 = ntRe(DiFfRG::zaaqbq2_qcd::tr15(fenv));
      const auto _ix9 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp62 = ZAqbq1.at(_ix9);
      const auto _interp64 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix17 = zq.index(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp65 = zq.at(_ix17);
      const auto _interp84 = ntRe(DiFfRG::zaaqbq2_qcd::tr17(fenv));
      const auto _interp85 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp86 = zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp87 = ntRe(DiFfRG::zaaqbq2_qcd::tr16(fenv));
      const auto _interp88 = ntRe(DiFfRG::zaaqbq2_qcd::tr18(fenv));
      const auto _ix1 = ZAqbq1.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp89 = ZAqbq1.at(_ix1);
      const auto _interp90 = ntRe(DiFfRG::zaaqbq2_qcd::tr19(fenv));
      const auto _ix2 = ZAqbq1.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp91 = ZAqbq1.at(_ix2);
      const auto _den14 = powr<-1>(_interp48 + _interp47 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den15 = powr<-1>(_interp58 + _interp57 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den16 = powr<-1>(_interp65 + _interp64 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den19 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den25 = 3. * powr<-1>(3. * _interp33 * nthk0 + _interp34 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den26 = 3. * powr<-1>(3. * _interp22 * nthk0 + _interp23 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den27 = powr<-1>(_interp50 * nthk0 + _interp51 * (powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den29 = powr<-1>(_interp86 + _interp85 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp1 = ntRe(DiFfRG::zaaqbq2_qcd::tr1(fenv));
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp6 = ZA.at(_ix13);
      const auto _interp9 = ZA3.at(_ix7);
      const auto _interp10 = ZA3.at(_ix1);
      const auto _interp14 = ntRe(DiFfRG::zaaqbq2_qcd::tr0(fenv));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp16 = ZA.at(_ix17);
      const auto _interp17 = ZA3.at(_ix9);
      const auto _interp18 = ZAAqbq2(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp66 = ntRe(DiFfRG::zaaqbq2_qcd::tr31(fenv));
      const auto _interp68 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp69 = zq.at(_ix11);
      const auto _interp70 = ntRe(DiFfRG::zaaqbq2_qcd::tr32(fenv));
      const auto _interp71 = ntRe(DiFfRG::zaaqbq2_qcd::tr33(fenv));
      const auto _interp72 = ZAAqbq2.at(_ix6);
      const auto _interp73 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp74 = zq.at(_ix18);
      const auto _interp75 = ntRe(DiFfRG::zaaqbq2_qcd::tr21(fenv));
      const auto _interp76 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp77 = ZA.at(_ix16);
      const auto _interp78 = ZA3.at(_ix8);
      const auto _interp79 = ZA3.at(_ix2);
      const auto _interp81 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp82 = zq.at(_ix19);
      const auto _interp83 = ntRe(DiFfRG::zaaqbq2_qcd::tr22(fenv));
      const auto _interp92 = ntRe(DiFfRG::zaaqbq2_qcd::tr34(fenv));
      const auto _interp93 = ZAqbq1.at(_ix3);
      const auto _interp94 = ZAqbq1.at(_ix0);
      const auto _interp95 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp96 = zq.at(_ix12);
      const auto _interp97 = ntRe(DiFfRG::zaaqbq2_qcd::tr23(fenv));
      const auto _interp98 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp99 = zq.at(_ix14);
      const auto _interp100 = ntRe(DiFfRG::zaaqbq2_qcd::tr24(fenv));
      const auto _interp101 = ntRe(DiFfRG::zaaqbq2_qcd::tr26(fenv));
      const auto _interp102 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp103 = ntRe(DiFfRG::zaaqbq2_qcd::tr35(fenv));
      const auto _interp104 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp105 = ntRe(DiFfRG::zaaqbq2_qcd::tr28(fenv));
      const auto _interp106 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp107 = zq.at(_ix15);
      const auto _den17 = powr<-1>(_interp74 + _interp73 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den18 = powr<-1>(_interp82 + _interp81 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den20 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den21 = powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den22 = 3. * powr<-1>(3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p));
      const auto _den30 = powr<-1>(_interp96 + _interp95 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den31 = powr<-1>(_interp99 + _interp98 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den32 = powr<-1>(_interp107 + _interp106 * sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den33 = powr<-1>(_interp69 + _interp68 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den20 + _den21 + _den22 + _den24 + _den25 + _den26 + _den27 + _den29 + _den3 + _den30 + _den31 + _den32 + _den33 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p1 + cosl1p2 + cosl1p3 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den2 = powr<-2>(_interp25 + _interp20 * powr<-1>(l1));
        const auto _den4 = powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den6 = powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den12 = powr<-1>(_interp29 * nthk0 + _interp30 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
        const auto _den13 = powr<-1>(_interp42 * nthk0 + _interp43 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
        const auto _den28 = powr<-1>(_interp31 * nthk0 + _interp32 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.000732421875, _den2 * _den26 * _interp19 * powr<2>(_interp24) * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.000732421875, _den2 * _den26 * _interp26 * powr<2>(_interp27) * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den12 * _den2 * _den25 * _den28 * _interp28 * _interp35 * _interp36 * _interp37 * _interp38 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.0003662109375, _den12 * _den2 * _den28 * _interp37 * _interp38 * _interp39 * _interp40 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.000732421875, _den13 * _den14 * _den2 * _den4 * _interp41 * _interp44 * _interp45 * _interp46 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.000732421875, _den14 * _den2 * _den27 * _den4 * _interp45 * _interp49 * _interp52 * _interp53 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den15 * _den2 * _den28 * _den5 * _interp38 * _interp54 * _interp55 * _interp56 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den15 * _den2 * _den26 * _den5 * _interp27 * _interp56 * _interp59 * _interp60 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den14 * _den16 * _den2 * _den27 * _den4 * _den6 * _interp45 * _interp53 * _interp61 * _interp62 * _interp63 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.0003662109375, _den19 * _den2 * _den28 * _den29 * _interp11 * _interp38 * _interp67 * _interp84 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.0003662109375, _den13 * _den19 * _den2 * _den29 * _interp11 * _interp46 * _interp80 * _interp87 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.000732421875, _den14 * _den19 * _den2 * _den28 * _den29 * _den4 * _interp38 * _interp45 * _interp67 * _interp88 * _interp89 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den13 * _den15 * _den19 * _den2 * _den29 * _den5 * _interp46 * _interp56 * _interp80 * _interp90 * _interp91 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.0003662109375, _den24 * _den3 * _interp11 * _interp12 * _interp13 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den7 = powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den8 = powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den9 = powr<-1>(_interp5 * nthk0 + _interp6 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp76 * nthk0 + _interp77 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den11 = powr<-1>(_interp15 * nthk0 + _interp16 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den23 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
        // clang-format off
        _acc += fma(-0.000732421875, _den21 * _den3 * _den31 * _interp100 * powr<2>(_interp27) * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den17 * _den21 * _den3 * _den31 * _den7 * _interp101 * _interp102 * _interp27 * _interp37 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den21 * _den23 * _den3 * _den31 * _den33 * _interp103 * _interp104 * _interp27 * _interp38 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den23 * _den24 * _den3 * _den33 * _interp13 * _interp38 * _interp67 * _interp70 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den17 * _den23 * _den3 * _den33 * _den7 * _interp37 * _interp38 * _interp71 * _interp72 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den10 * _den18 * _den24 * _den3 * _den8 * _interp46 * _interp75 * _interp78 * _interp79 * _interp80 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den18 * _den24 * _den3 * _den8 * _interp13 * _interp46 * _interp80 * _interp83 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den24 * _den3 * _den9 * _interp1 * _interp10 * _interp11 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.0003662109375, _den11 * _den3 * _den9 * _interp14 * _interp17 * _interp18 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den11 * _den22 * _den3 * _den32 * _den9 * _interp105 * _interp17 * _interp53 * _interp63 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den23 * _den24 * _den3 * _den33 * _den9 * _interp10 * _interp38 * _interp66 * _interp67 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den17 * _den20 * _den23 * _den3 * _den30 * _den33 * _den7 * _interp37 * _interp38 * _interp92 * _interp93 * _interp94 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den21 * _den3 * _den31 * powr<2>(_interp24) * _interp97 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), 0.)))))))))))));
        // clang-format on
      }
      return _acc;
    }
    #elif NT_ZAAQBQ2_QCD_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaaqbq2_qcd::nenv) > 0 ? (DiFfRG::zaaqbq2_qcd::nenv) : 1];
      DiFfRG::zaaqbq2_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix12 = ZA.index(l1);
      const auto _interp4 = ZA.at(_ix12);
      const auto _interp7 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp8 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _ix7 = ZAAqbq2.index(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp11 = ZAAqbq2.at(_ix7);
      const auto _interp13 = ZA4.at(_ix7);
      const auto _interp20 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp21 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _ix8 = ZAAqbq1.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp24 = ZAAqbq1.at(_ix8);
      const auto _interp25 = zq.at(_ix12);
      const auto _interp27 = ZAAqbq2.at(_ix8);
      const auto _interp32 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp34 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _ix13 = ZA.index(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp35 = ZA.at(_ix13);
      const auto _interp37 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp46 = RB(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp47 = ZA(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp48 = ZAAqbq2(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p3, l1 * p, fma(2., powr<2>(p), powr<2>(l1))))));
      const auto _interp49 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp50 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp51 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp52 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp56 = ZAAqbq2(0.7071067811865475 * sqrt(fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp57 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp61 = ZAAqbq2(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(2., powr<2>(p), powr<2>(l1)))))));
      const auto _ix10 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp62 = ZAqbq1.at(_ix10);
      const auto _interp63 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _ix14 = zq.index(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp64 = zq.at(_ix14);
      const auto _interp66 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp70 = ZAqbq1(0.816496580927726 * sqrt(fma(2., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp74 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp75 = ZA.at(_ix14);
      const auto _interp76 = ZA3.at(_ix10);
      const auto _interp77 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp78 = zq.at(_ix13);
      const auto _interp80 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp84 = RF(powr<2>(k), fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp85 = zq(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp91 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp20, powr<-1>(l1), _interp25));
      const auto _den3 = powr<-2>(fma(_interp4, powr<2>(l1), fma(_interp2, nthk0, 0.)));
      const auto _den4 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den5 = powr<-1>(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den7 = powr<-1>(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den10 = powr<-1>(fma(_interp74, nthk0, fma(_interp75, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den13 = powr<-1>(fma(_interp46, nthk0, fma(_interp47, powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den14 = powr<-1>(fma(_interp51, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp52));
      const auto _den15 = powr<-1>(fma(_interp63, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))), _interp64));
      const auto _den17 = powr<-1>(fma(_interp84, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))), _interp85));
      const auto _den23 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _den24 = powr<-1>(fma(_interp7, nthk0, fma(_interp8, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den28 = powr<-1>(fma(_interp34, nthk0, fma(_interp35, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
      const auto _den33 = powr<-1>(fma(_interp77, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))), _interp78));
      const auto _interp19 = ntRe(DiFfRG::zaaqbq2_qcd::tr3(fenv));
      const auto _interp22 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix16 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA.at(_ix16);
      const auto _interp26 = ntRe(DiFfRG::zaaqbq2_qcd::tr4(fenv));
      const auto _interp28 = ntIm(DiFfRG::zaaqbq2_qcd::tr5(fenv));
      const auto _interp29 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp30 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _ix5 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp31 = ZA3.at(_ix5);
      const auto _interp33 = ntIm(DiFfRG::zaaqbq2_qcd::tr6(fenv));
      const auto _ix2 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp36 = ZA3.at(_ix2);
      const auto _interp38 = ntRe(DiFfRG::zaaqbq2_qcd::tr7(fenv));
      const auto _interp39 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix15 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp40 = ZA.at(_ix15);
      const auto _ix4 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp41 = ZA3.at(_ix4);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp42 = ZA3.at(_ix0);
      const auto _interp43 = ntRe(DiFfRG::zaaqbq2_qcd::tr8(fenv));
      const auto _ix9 = ZA4.index(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp44 = ZA4.at(_ix9);
      const auto _interp45 = ntRe(DiFfRG::zaaqbq2_qcd::tr9(fenv));
      const auto _interp53 = ntRe(DiFfRG::zaaqbq2_qcd::tr10(fenv));
      const auto _interp54 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix17 = ZA.index(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp55 = ZA.at(_ix17);
      const auto _interp58 = ntIm(DiFfRG::zaaqbq2_qcd::tr11(fenv));
      const auto _ix6 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp59 = ZA3.at(_ix6);
      const auto _interp60 = ntRe(DiFfRG::zaaqbq2_qcd::tr13(fenv));
      const auto _interp65 = ntRe(DiFfRG::zaaqbq2_qcd::tr12(fenv));
      const auto _interp67 = ntIm(DiFfRG::zaaqbq2_qcd::tr14(fenv));
      const auto _interp68 = ntRe(DiFfRG::zaaqbq2_qcd::tr15(fenv));
      const auto _ix11 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp69 = ZAqbq1.at(_ix11);
      const auto _interp71 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix18 = zq.index(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp72 = zq.at(_ix18);
      const auto _den16 = powr<-1>(_interp72 + _interp71 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den25 = 3. * powr<-1>(3. * _interp39 * nthk0 + _interp40 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den26 = 3. * powr<-1>(3. * _interp22 * nthk0 + _interp23 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den27 = powr<-1>(_interp54 * nthk0 + _interp55 * (powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp12 = ntRe(DiFfRG::zaaqbq2_qcd::tr2(fenv));
      const auto _interp73 = ntIm(DiFfRG::zaaqbq2_qcd::tr30(fenv));
      const auto _interp93 = ntRe(DiFfRG::zaaqbq2_qcd::tr17(fenv));
      const auto _interp94 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp95 = zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp96 = ntRe(DiFfRG::zaaqbq2_qcd::tr16(fenv));
      const auto _interp97 = ntRe(DiFfRG::zaaqbq2_qcd::tr18(fenv));
      const auto _ix1 = ZAqbq1.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp98 = ZAqbq1.at(_ix1);
      const auto _interp99 = ntRe(DiFfRG::zaaqbq2_qcd::tr19(fenv));
      const auto _ix3 = ZAqbq1.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp100 = ZAqbq1.at(_ix3);
      const auto _interp101 = ntRe(DiFfRG::zaaqbq2_qcd::tr34(fenv));
      const auto _interp102 = ZAqbq1.at(_ix4);
      const auto _interp103 = ZAqbq1.at(_ix0);
      const auto _interp104 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp105 = zq.at(_ix15);
      const auto _interp106 = ntRe(DiFfRG::zaaqbq2_qcd::tr23(fenv));
      const auto _interp107 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp108 = zq.at(_ix16);
      const auto _interp109 = ntRe(DiFfRG::zaaqbq2_qcd::tr24(fenv));
      const auto _interp110 = ntIm(DiFfRG::zaaqbq2_qcd::tr25(fenv));
      const auto _interp111 = ntRe(DiFfRG::zaaqbq2_qcd::tr26(fenv));
      const auto _interp112 = ZAqbq1.at(_ix5);
      const auto _interp113 = ntRe(DiFfRG::zaaqbq2_qcd::tr35(fenv));
      const auto _interp114 = ZAqbq1.at(_ix2);
      const auto _interp115 = ntIm(DiFfRG::zaaqbq2_qcd::tr36(fenv));
      const auto _den19 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den20 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den21 = powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den29 = powr<-1>(_interp95 + _interp94 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den30 = powr<-1>(_interp105 + _interp104 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den31 = powr<-1>(_interp108 + _interp107 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp1 = ntRe(DiFfRG::zaaqbq2_qcd::tr1(fenv));
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp6 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3.at(_ix1);
      const auto _interp14 = ntRe(DiFfRG::zaaqbq2_qcd::tr0(fenv));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp16 = ZA.at(_ix18);
      const auto _interp17 = ZA3.at(_ix11);
      const auto _interp18 = ZAAqbq2(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp79 = ntRe(DiFfRG::zaaqbq2_qcd::tr31(fenv));
      const auto _interp81 = ntRe(DiFfRG::zaaqbq2_qcd::tr32(fenv));
      const auto _interp82 = ntRe(DiFfRG::zaaqbq2_qcd::tr33(fenv));
      const auto _interp83 = ZAAqbq2.at(_ix9);
      const auto _interp86 = ntIm(DiFfRG::zaaqbq2_qcd::tr20(fenv));
      const auto _interp87 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp88 = zq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp89 = ntRe(DiFfRG::zaaqbq2_qcd::tr21(fenv));
      const auto _interp90 = ZA3.at(_ix3);
      const auto _interp92 = ntRe(DiFfRG::zaaqbq2_qcd::tr22(fenv));
      const auto _interp116 = ntIm(DiFfRG::zaaqbq2_qcd::tr27(fenv));
      const auto _interp117 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp118 = zq.at(_ix17);
      const auto _interp119 = ntRe(DiFfRG::zaaqbq2_qcd::tr28(fenv));
      const auto _interp120 = ntIm(DiFfRG::zaaqbq2_qcd::tr29(fenv));
      const auto _interp121 = ZAqbq1.at(_ix6);
      const auto _den18 = powr<-1>(_interp88 + _interp87 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den22 = 3. * powr<-1>(3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p));
      const auto _den32 = powr<-1>(_interp118 + _interp117 * sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den2 + _den20 + _den21 + _den22 + _den23 + _den24 + _den25 + _den26 + _den27 + _den28 + _den29 + _den3 + _den30 + _den31 + _den32 + _den33 + _den4 + _den5 + _den7 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p1 + cosl1p2 + cosl1p3 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den6 = powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den12 = powr<-1>(_interp29 * nthk0 + _interp30 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += -0.0003662109375 * fma(-2., _den2 * _den26 * _interp19 * _interp21 * powr<2>(_interp24) * powr<-5>(l1) * powr<-2>(p), fma(2., _den2 * _den26 * _interp21 * _interp26 * powr<2>(_interp27) * powr<-5>(l1) * powr<-2>(p), fma(2., _den12 * _den2 * _den26 * _interp21 * _interp27 * _interp28 * _interp31 * _interp32 * powr<-5>(l1) * powr<-2>(p), fma(2., _den2 * _den26 * _den28 * _interp21 * _interp27 * _interp33 * _interp36 * _interp37 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den12 * _den2 * _den25 * _den28 * _interp21 * _interp32 * _interp37 * _interp38 * _interp41 * _interp42 * powr<-5>(l1) * powr<-2>(p), fma(-1., _den12 * _den2 * _den28 * _interp21 * _interp32 * _interp37 * _interp43 * _interp44 * powr<-5>(l1) * powr<-2>(p), fma(2., _den13 * _den14 * _den2 * _den4 * _interp21 * _interp45 * _interp48 * _interp49 * _interp50 * powr<-5>(l1) * powr<-2>(p), fma(2., _den14 * _den2 * _den27 * _den4 * _interp21 * _interp49 * _interp53 * _interp56 * _interp57 * powr<-5>(l1) * powr<-2>(p), fma(2., _den13 * _den14 * _den2 * _den27 * _den4 * _interp21 * _interp49 * _interp50 * _interp57 * _interp58 * _interp59 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den15 * _den2 * _den28 * _den5 * _interp21 * _interp37 * _interp60 * _interp61 * _interp62 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den15 * _den2 * _den26 * _den5 * _interp21 * _interp27 * _interp62 * _interp65 * _interp66 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den15 * _den2 * _den26 * _den28 * _den5 * _interp21 * _interp36 * _interp37 * _interp62 * _interp66 * _interp67 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den14 * _den16 * _den2 * _den27 * _den4 * _den6 * _interp21 * _interp49 * _interp57 * _interp68 * _interp69 * _interp70 * powr<-5>(l1) * powr<-2>(p), fma(2., _den2 * _den26 * _interp19 * _interp20 * powr<2>(_interp24) * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den2 * _den26 * _interp20 * _interp26 * powr<2>(_interp27) * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den12 * _den2 * _den26 * _interp20 * _interp27 * _interp28 * _interp31 * _interp32 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den2 * _den26 * _den28 * _interp20 * _interp27 * _interp33 * _interp36 * _interp37 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den12 * _den2 * _den25 * _den28 * _interp20 * _interp32 * _interp37 * _interp38 * _interp41 * _interp42 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(_den12, _den2 * _den28 * _interp20 * _interp32 * _interp37 * _interp43 * _interp44 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den13 * _den14 * _den2 * _den4 * _interp20 * _interp45 * _interp48 * _interp49 * _interp50 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den14 * _den2 * _den27 * _den4 * _interp20 * _interp49 * _interp53 * _interp56 * _interp57 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den13 * _den14 * _den2 * _den27 * _den4 * _interp20 * _interp49 * _interp50 * _interp57 * _interp58 * _interp59 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den15 * _den2 * _den28 * _den5 * _interp20 * _interp37 * _interp60 * _interp61 * _interp62 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den15 * _den2 * _den26 * _den5 * _interp20 * _interp27 * _interp62 * _interp65 * _interp66 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den15 * _den2 * _den26 * _den28 * _den5 * _interp20 * _interp36 * _interp37 * _interp62 * _interp66 * _interp67 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den14 * _den16 * _den2 * _den27 * _den4 * _den6 * _interp20 * _interp49 * _interp57 * _interp68 * _interp69 * _interp70 * etaQ * powr<-5>(l1) * powr<-2>(p), 0.))))))))))))))))))))))))));
        // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(0.0003662109375, _den19 * _den2 * _den28 * _den29 * _interp11 * _interp37 * _interp80 * _interp93 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.0003662109375, _den13 * _den19 * _den2 * _den29 * _interp11 * _interp50 * _interp91 * _interp96 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.000732421875, _den14 * _den19 * _den2 * _den28 * _den29 * _den4 * _interp37 * _interp49 * _interp80 * _interp97 * _interp98 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den13 * _den15 * _den19 * _den2 * _den29 * _den5 * _interp100 * _interp50 * _interp62 * _interp91 * _interp99 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.0003662109375, _den24 * _den3 * _interp11 * _interp12 * _interp13 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den21 * _den3 * _den31 * _interp106 * powr<2>(_interp24) * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den21 * _den3 * _den31 * _interp109 * powr<2>(_interp27) * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den17 * _den21 * _den3 * _den31 * _den7 * _interp111 * _interp112 * _interp27 * _interp32 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den21 * _den23 * _den3 * _den31 * _den33 * _interp113 * _interp114 * _interp27 * _interp37 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den17 * _den20 * _den23 * _den3 * _den30 * _den33 * _den7 * _interp101 * _interp102 * _interp103 * _interp32 * _interp37 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den10 * _den21 * _den3 * _den31 * _interp110 * _interp27 * _interp66 * _interp76 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den10 * _den21 * _den23 * _den3 * _den31 * _den33 * _interp114 * _interp115 * _interp37 * _interp66 * _interp76 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den10 * _den23 * _den3 * _den33 * _interp37 * _interp61 * _interp73 * _interp76 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), 0.)))))))))))));
      // clang-format on
      }
      { // subkernel 3
        const auto _den8 = powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den9 = powr<-1>(_interp5 * nthk0 + _interp6 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den11 = powr<-1>(_interp15 * nthk0 + _interp16 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.0003662109375, _den23 * _den24 * _den3 * _den33 * _interp13 * _interp37 * _interp80 * _interp81 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den17 * _den23 * _den3 * _den33 * _den7 * _interp32 * _interp37 * _interp82 * _interp83 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den24 * _den3 * _den9 * _interp1 * _interp10 * _interp11 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.0003662109375, _den11 * _den3 * _den9 * _interp14 * _interp17 * _interp18 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den18 * _den22 * _den3 * _den32 * _den8 * _den9 * _interp120 * _interp121 * _interp50 * _interp57 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den22 * _den3 * _den32 * _den9 * _interp116 * _interp56 * _interp57 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den11 * _den22 * _den3 * _den32 * _den9 * _interp119 * _interp17 * _interp57 * _interp70 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den23 * _den24 * _den3 * _den33 * _den9 * _interp10 * _interp37 * _interp79 * _interp80 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den18 * _den3 * _den8 * _den9 * _interp48 * _interp50 * _interp86 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den10 * _den18 * _den24 * _den3 * _den8 * _interp50 * _interp76 * _interp89 * _interp90 * _interp91 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den18 * _den24 * _den3 * _den8 * _interp13 * _interp50 * _interp91 * _interp92 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), 0.)))))))))));
        // clang-format on
      }
      return _acc;
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaaqbq2_qcd::nenv) > 0 ? (DiFfRG::zaaqbq2_qcd::nenv) : 1];
      DiFfRG::zaaqbq2_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix11 = ZA.index(l1);
      const auto _interp3 = ZA.at(_ix11);
      const auto _interp6 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp7 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp10 = ZAAqbq2(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp16 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp17 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _ix7 = ZAAqbq1.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp20 = ZAAqbq1.at(_ix7);
      const auto _interp21 = zq.at(_ix11);
      const auto _interp22 = ZAAqbq2.at(_ix7);
      const auto _interp26 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _ix12 = ZA.index(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp28 = ZA.at(_ix12);
      const auto _interp30 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp36 = RB(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _ix14 = ZA.index(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp37 = ZA.at(_ix14);
      const auto _interp38 = ZAAqbq2(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p3, l1 * p, fma(2., powr<2>(p), powr<2>(l1))))));
      const auto _interp39 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp40 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp41 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp42 = zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp45 = ZAAqbq2(0.7071067811865475 * sqrt(fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp46 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp48 = ZAAqbq2(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(2., powr<2>(p), powr<2>(l1)))))));
      const auto _ix9 = ZAqbq1.index(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp49 = ZAqbq1.at(_ix9);
      const auto _interp50 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _ix13 = zq.index(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp51 = zq.at(_ix13);
      const auto _interp52 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp54 = ZAqbq1(0.816496580927726 * sqrt(fma(2., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp57 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp58 = ZA.at(_ix13);
      const auto _interp59 = ZA3.at(_ix9);
      const auto _interp60 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp61 = zq.at(_ix12);
      const auto _interp62 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp64 = RF(powr<2>(k), fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp65 = zq(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp66 = RF(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp67 = zq.at(_ix14);
      const auto _interp78 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp79 = zq(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp16, powr<-1>(l1), _interp21));
      const auto _den3 = powr<-2>(fma(_interp3, powr<2>(l1), fma(_interp1, nthk0, 0.)));
      const auto _den4 = powr<-1>(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den5 = powr<-1>(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den7 = powr<-1>(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den8 = powr<-1>(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _den10 = powr<-1>(fma(_interp57, nthk0, fma(_interp58, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den13 = powr<-1>(fma(_interp36, nthk0, fma(_interp37, powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den14 = powr<-1>(fma(_interp41, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))), _interp42));
      const auto _den15 = powr<-1>(fma(_interp50, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))), _interp51));
      const auto _den17 = powr<-1>(fma(_interp64, sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))), _interp65));
      const auto _den18 = powr<-1>(fma(_interp66, sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))), _interp67));
      const auto _den21 = powr<-1>(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _den23 = powr<-1>(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _den24 = powr<-1>(fma(_interp6, nthk0, fma(_interp7, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den28 = powr<-1>(fma(_interp27, nthk0, fma(_interp28, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
      const auto _den31 = powr<-1>(fma(_interp78, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))), _interp79));
      const auto _den33 = powr<-1>(fma(_interp60, sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))), _interp61));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp24 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _ix5 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp25 = ZA3.at(_ix5);
      const auto _ix2 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp29 = ZA3.at(_ix2);
      const auto _interp31 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix15 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp32 = ZA.at(_ix15);
      const auto _ix4 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp33 = ZA3.at(_ix4);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp34 = ZA3.at(_ix0);
      const auto _ix8 = ZA4.index(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp35 = ZA4.at(_ix8);
      const auto _interp43 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix16 = ZA.index(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp44 = ZA.at(_ix16);
      const auto _ix6 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp47 = ZA3.at(_ix6);
      const auto _ix10 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp53 = ZAqbq1.at(_ix10);
      const auto _interp55 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix17 = zq.index(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp56 = zq.at(_ix17);
      const auto _den16 = powr<-1>(_interp56 + _interp55 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den25 = 3. * powr<-1>(3. * _interp31 * nthk0 + _interp32 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den26 = 3. * powr<-1>(3. * _interp18 * nthk0 + _interp19 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den27 = powr<-1>(_interp43 * nthk0 + _interp44 * (powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp11 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp63 = ZAAqbq2.at(_ix8);
      const auto _ix3 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp68 = ZA3.at(_ix3);
      const auto _interp69 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp70 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp71 = zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _ix1 = ZAqbq1.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp72 = ZAqbq1.at(_ix1);
      const auto _interp73 = ZAqbq1.at(_ix3);
      const auto _den19 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den29 = powr<-1>(_interp71 + _interp70 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp5 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp8 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3.at(_ix1);
      const auto _interp12 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp13 = ZA.at(_ix17);
      const auto _interp14 = ZA3.at(_ix10);
      const auto _interp15 = ZAAqbq2(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp74 = ZAqbq1.at(_ix4);
      const auto _interp75 = ZAqbq1.at(_ix0);
      const auto _interp76 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp77 = zq.at(_ix15);
      const auto _interp80 = ZAqbq1.at(_ix5);
      const auto _interp81 = ZAqbq1.at(_ix2);
      const auto _interp82 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp83 = zq.at(_ix16);
      const auto _interp84 = ZAqbq1.at(_ix6);
      const auto _den20 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den22 = 3. * powr<-1>(3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p));
      const auto _den30 = powr<-1>(_interp77 + _interp76 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den32 = powr<-1>(_interp83 + _interp82 * sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      // clang-format off
      using _T = decltype(complex<double>(0.,-2.) + _den1 + _den10 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den2 + _den20 + _den21 + _den22 + _den23 + _den24 + _den25 + _den26 + _den27 + _den28 + _den29 + _den3 + _den30 + _den31 + _den32 + _den33 + _den4 + _den5 + _den7 + _den8 + DiFfRG::zaaqbq2_qcd::tr0(fenv) + DiFfRG::zaaqbq2_qcd::tr10(fenv) + DiFfRG::zaaqbq2_qcd::tr11(fenv) + DiFfRG::zaaqbq2_qcd::tr12(fenv) + DiFfRG::zaaqbq2_qcd::tr13(fenv) + DiFfRG::zaaqbq2_qcd::tr14(fenv) + DiFfRG::zaaqbq2_qcd::tr15(fenv) + DiFfRG::zaaqbq2_qcd::tr16(fenv) + DiFfRG::zaaqbq2_qcd::tr17(fenv) + DiFfRG::zaaqbq2_qcd::tr18(fenv) + DiFfRG::zaaqbq2_qcd::tr19(fenv) + DiFfRG::zaaqbq2_qcd::tr1(fenv) + DiFfRG::zaaqbq2_qcd::tr20(fenv) + DiFfRG::zaaqbq2_qcd::tr21(fenv) + DiFfRG::zaaqbq2_qcd::tr22(fenv) + DiFfRG::zaaqbq2_qcd::tr23(fenv) + DiFfRG::zaaqbq2_qcd::tr24(fenv) + DiFfRG::zaaqbq2_qcd::tr25(fenv) + DiFfRG::zaaqbq2_qcd::tr26(fenv) + DiFfRG::zaaqbq2_qcd::tr27(fenv) + DiFfRG::zaaqbq2_qcd::tr28(fenv) + DiFfRG::zaaqbq2_qcd::tr29(fenv) + DiFfRG::zaaqbq2_qcd::tr2(fenv) + DiFfRG::zaaqbq2_qcd::tr30(fenv) + DiFfRG::zaaqbq2_qcd::tr31(fenv) + DiFfRG::zaaqbq2_qcd::tr32(fenv) + DiFfRG::zaaqbq2_qcd::tr33(fenv) + DiFfRG::zaaqbq2_qcd::tr34(fenv) + DiFfRG::zaaqbq2_qcd::tr35(fenv) + DiFfRG::zaaqbq2_qcd::tr36(fenv) + DiFfRG::zaaqbq2_qcd::tr3(fenv) + DiFfRG::zaaqbq2_qcd::tr4(fenv) + DiFfRG::zaaqbq2_qcd::tr5(fenv) + DiFfRG::zaaqbq2_qcd::tr6(fenv) + DiFfRG::zaaqbq2_qcd::tr7(fenv) + DiFfRG::zaaqbq2_qcd::tr8(fenv) + DiFfRG::zaaqbq2_qcd::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den6 = powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den12 = powr<-1>(_interp23 * nthk0 + _interp24 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += -0.0003662109375 * fma(-2., _den2 * _den26 * DiFfRG::zaaqbq2_qcd::tr3(fenv) * _interp17 * powr<2>(_interp20) * powr<-5>(l1) * powr<-2>(p), fma(2., _den2 * _den26 * DiFfRG::zaaqbq2_qcd::tr4(fenv) * _interp17 * powr<2>(_interp22) * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,-2.), _den12 * _den2 * _den26 * DiFfRG::zaaqbq2_qcd::tr5(fenv) * _interp17 * _interp22 * _interp25 * _interp26 * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,-2.), _den2 * _den26 * _den28 * DiFfRG::zaaqbq2_qcd::tr6(fenv) * _interp17 * _interp22 * _interp29 * _interp30 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den12 * _den2 * _den25 * _den28 * DiFfRG::zaaqbq2_qcd::tr7(fenv) * _interp17 * _interp26 * _interp30 * _interp33 * _interp34 * powr<-5>(l1) * powr<-2>(p), fma(-1., _den12 * _den2 * _den28 * DiFfRG::zaaqbq2_qcd::tr8(fenv) * _interp17 * _interp26 * _interp30 * _interp35 * powr<-5>(l1) * powr<-2>(p), fma(2., _den13 * _den14 * _den2 * _den4 * DiFfRG::zaaqbq2_qcd::tr9(fenv) * _interp17 * _interp38 * _interp39 * _interp40 * powr<-5>(l1) * powr<-2>(p), fma(2., _den14 * _den2 * _den27 * _den4 * DiFfRG::zaaqbq2_qcd::tr10(fenv) * _interp17 * _interp39 * _interp45 * _interp46 * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,-2.), _den13 * _den14 * _den2 * _den27 * _den4 * DiFfRG::zaaqbq2_qcd::tr11(fenv) * _interp17 * _interp39 * _interp40 * _interp46 * _interp47 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den15 * _den2 * _den28 * _den5 * DiFfRG::zaaqbq2_qcd::tr13(fenv) * _interp17 * _interp30 * _interp48 * _interp49 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den15 * _den2 * _den26 * _den5 * DiFfRG::zaaqbq2_qcd::tr12(fenv) * _interp17 * _interp22 * _interp49 * _interp52 * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,2.), _den15 * _den2 * _den26 * _den28 * _den5 * DiFfRG::zaaqbq2_qcd::tr14(fenv) * _interp17 * _interp29 * _interp30 * _interp49 * _interp52 * powr<-5>(l1) * powr<-2>(p), fma(-2., _den14 * _den16 * _den2 * _den27 * _den4 * _den6 * DiFfRG::zaaqbq2_qcd::tr15(fenv) * _interp17 * _interp39 * _interp46 * _interp53 * _interp54 * powr<-5>(l1) * powr<-2>(p), fma(2., _den2 * _den26 * DiFfRG::zaaqbq2_qcd::tr3(fenv) * _interp16 * powr<2>(_interp20) * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den2 * _den26 * DiFfRG::zaaqbq2_qcd::tr4(fenv) * _interp16 * powr<2>(_interp22) * etaQ * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,2.), _den12 * _den2 * _den26 * DiFfRG::zaaqbq2_qcd::tr5(fenv) * _interp16 * _interp22 * _interp25 * _interp26 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,2.), _den2 * _den26 * _den28 * DiFfRG::zaaqbq2_qcd::tr6(fenv) * _interp16 * _interp22 * _interp29 * _interp30 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den12 * _den2 * _den25 * _den28 * DiFfRG::zaaqbq2_qcd::tr7(fenv) * _interp16 * _interp26 * _interp30 * _interp33 * _interp34 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(_den12, _den2 * _den28 * DiFfRG::zaaqbq2_qcd::tr8(fenv) * _interp16 * _interp26 * _interp30 * _interp35 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den13 * _den14 * _den2 * _den4 * DiFfRG::zaaqbq2_qcd::tr9(fenv) * _interp16 * _interp38 * _interp39 * _interp40 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(-2., _den14 * _den2 * _den27 * _den4 * DiFfRG::zaaqbq2_qcd::tr10(fenv) * _interp16 * _interp39 * _interp45 * _interp46 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,2.), _den13 * _den14 * _den2 * _den27 * _den4 * DiFfRG::zaaqbq2_qcd::tr11(fenv) * _interp16 * _interp39 * _interp40 * _interp46 * _interp47 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den15 * _den2 * _den28 * _den5 * DiFfRG::zaaqbq2_qcd::tr13(fenv) * _interp16 * _interp30 * _interp48 * _interp49 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den15 * _den2 * _den26 * _den5 * DiFfRG::zaaqbq2_qcd::tr12(fenv) * _interp16 * _interp22 * _interp49 * _interp52 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(complex<double>(0.,-2.), _den15 * _den2 * _den26 * _den28 * _den5 * DiFfRG::zaaqbq2_qcd::tr14(fenv) * _interp16 * _interp29 * _interp30 * _interp49 * _interp52 * etaQ * powr<-5>(l1) * powr<-2>(p), fma(2., _den14 * _den16 * _den2 * _den27 * _den4 * _den6 * DiFfRG::zaaqbq2_qcd::tr15(fenv) * _interp16 * _interp39 * _interp46 * _interp53 * _interp54 * etaQ * powr<-5>(l1) * powr<-2>(p), 0.))))))))))))))))))))))))));
        // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(0.0003662109375, _den19 * _den2 * _den28 * _den29 * DiFfRG::zaaqbq2_qcd::tr17(fenv) * _interp10 * _interp30 * _interp62 * (_interp17 - _interp16 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.0003662109375, _den13 * _den19 * _den2 * _den29 * DiFfRG::zaaqbq2_qcd::tr16(fenv) * _interp10 * _interp40 * _interp69 * (_interp17 - _interp16 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.000732421875, _den14 * _den19 * _den2 * _den28 * _den29 * _den4 * DiFfRG::zaaqbq2_qcd::tr18(fenv) * _interp30 * _interp39 * _interp62 * _interp72 * (_interp17 - _interp16 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(0.000732421875, _den13 * _den15 * _den19 * _den2 * _den29 * _den5 * DiFfRG::zaaqbq2_qcd::tr19(fenv) * _interp40 * _interp49 * _interp69 * _interp73 * (_interp17 - _interp16 * etaQ) * powr<-5>(l1) * powr<-2>(p), fma(-0.0003662109375, _den24 * _den3 * DiFfRG::zaaqbq2_qcd::tr2(fenv) * _interp10 * _interp11 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den21 * _den3 * _den31 * DiFfRG::zaaqbq2_qcd::tr23(fenv) * powr<2>(_interp20) * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den21 * _den3 * _den31 * DiFfRG::zaaqbq2_qcd::tr24(fenv) * powr<2>(_interp22) * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(complex<double>(0.,0.000732421875), _den10 * _den23 * _den3 * _den33 * DiFfRG::zaaqbq2_qcd::tr30(fenv) * _interp30 * _interp48 * _interp59 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(complex<double>(0.,0.000732421875), _den10 * _den21 * _den3 * _den31 * DiFfRG::zaaqbq2_qcd::tr25(fenv) * _interp22 * _interp52 * _interp59 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den23 * _den24 * _den3 * _den33 * DiFfRG::zaaqbq2_qcd::tr32(fenv) * _interp11 * _interp30 * _interp62 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den17 * _den23 * _den3 * _den33 * _den7 * DiFfRG::zaaqbq2_qcd::tr33(fenv) * _interp26 * _interp30 * _interp63 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.0003662109375, _den18 * _den24 * _den3 * _den8 * DiFfRG::zaaqbq2_qcd::tr22(fenv) * _interp11 * _interp40 * _interp69 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den10 * _den18 * _den24 * _den3 * _den8 * DiFfRG::zaaqbq2_qcd::tr21(fenv) * _interp40 * _interp59 * _interp68 * _interp69 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), 0.)))))))))))));
      // clang-format on
      }
      { // subkernel 3
        const auto _den9 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den11 = powr<-1>(_interp12 * nthk0 + _interp13 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(-0.000732421875, _den17 * _den20 * _den23 * _den3 * _den30 * _den33 * _den7 * DiFfRG::zaaqbq2_qcd::tr34(fenv) * _interp26 * _interp30 * _interp74 * _interp75 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.0003662109375, _den11 * _den3 * _den9 * DiFfRG::zaaqbq2_qcd::tr0(fenv) * _interp14 * _interp15 * _interp8 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(complex<double>(0.,0.000732421875), _den18 * _den3 * _den8 * _den9 * DiFfRG::zaaqbq2_qcd::tr20(fenv) * _interp38 * _interp40 * _interp8 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(complex<double>(0.,0.000732421875), _den22 * _den3 * _den32 * _den9 * DiFfRG::zaaqbq2_qcd::tr27(fenv) * _interp45 * _interp46 * _interp8 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den11 * _den22 * _den3 * _den32 * _den9 * DiFfRG::zaaqbq2_qcd::tr28(fenv) * _interp14 * _interp46 * _interp54 * _interp8 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den17 * _den21 * _den3 * _den31 * _den7 * DiFfRG::zaaqbq2_qcd::tr26(fenv) * _interp22 * _interp26 * _interp80 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den21 * _den23 * _den3 * _den31 * _den33 * DiFfRG::zaaqbq2_qcd::tr35(fenv) * _interp22 * _interp30 * _interp81 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(complex<double>(0.,-0.000732421875), _den10 * _den21 * _den23 * _den3 * _den31 * _den33 * DiFfRG::zaaqbq2_qcd::tr36(fenv) * _interp30 * _interp52 * _interp59 * _interp81 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(complex<double>(0.,-0.000732421875), _den18 * _den22 * _den3 * _den32 * _den8 * _den9 * DiFfRG::zaaqbq2_qcd::tr29(fenv) * _interp40 * _interp46 * _interp8 * _interp84 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(-0.000732421875, _den24 * _den3 * _den9 * DiFfRG::zaaqbq2_qcd::tr1(fenv) * _interp10 * _interp8 * _interp9 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), fma(0.000732421875, _den23 * _den24 * _den3 * _den33 * _den9 * DiFfRG::zaaqbq2_qcd::tr31(fenv) * _interp30 * _interp62 * _interp8 * _interp9 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))) * powr<-2>(p), 0.)))))))))));
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
      const auto _interp1 = ZAAqbq2(p);
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
using DiFfRG::ZAAqbq2_kernel;