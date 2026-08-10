#pragma once

#include "DiFfRG/physics/interpolation.hh"
#include "DiFfRG/physics/physics.hh"
#include "numtracer/sun/sun_data.hpp"
#include "kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  template<typename REG>
  class ZAAqbq1_kernel
  {
    public:
    using Regulator = REG;

    #if NT_ZAAQBQ1_QCD_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaaqbq1_qcd::nenv) > 0 ? (DiFfRG::zaaqbq1_qcd::nenv) : 1];
      DiFfRG::zaaqbq1_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp1 = ntRe(DiFfRG::zaaqbq1_qcd::tr1(fenv));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix3 = ZA.index(l1);
      const auto _interp4 = ZA.at(_ix3);
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix6 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp6 = ZA.at(_ix6);
      const auto _interp7 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix5 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp8 = ZA.at(_ix5);
      const auto _ix2 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp9 = ZA3.at(_ix2);
      const auto _interp10 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _ix0 = ZAAqbq1.index(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp11 = ZAAqbq1.at(_ix0);
      const auto _interp12 = ntRe(DiFfRG::zaaqbq1_qcd::tr2(fenv));
      const auto _interp13 = ZA4.at(_ix0);
      const auto _interp14 = ntRe(DiFfRG::zaaqbq1_qcd::tr0(fenv));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp16 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp18 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp19 = ntRe(DiFfRG::zaaqbq1_qcd::tr3(fenv));
      const auto _interp20 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp21 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp22 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix7 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA.at(_ix7);
      const auto _ix1 = ZAAqbq1.index(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp24 = ZAAqbq1.at(_ix1);
      const auto _interp25 = ZAAqbq2.at(_ix1);
      const auto _interp26 = zq.at(_ix3);
      const auto _interp27 = ntRe(DiFfRG::zaaqbq1_qcd::tr6(fenv));
      const auto _interp28 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp29 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp30 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp31 = ZAqbq1.at(_ix2);
      const auto _interp32 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp33 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp34 = zq.at(_ix6);
      const auto _interp35 = ntRe(DiFfRG::zaaqbq1_qcd::tr7(fenv));
      const auto _interp36 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp37 = ZA(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp38 = ZAAqbq1(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp39 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp40 = ntRe(DiFfRG::zaaqbq1_qcd::tr9(fenv));
      const auto _interp41 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _ix4 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp42 = ZA.at(_ix4);
      const auto _interp43 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp44 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp45 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp46 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp47 = zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp48 = ntRe(DiFfRG::zaaqbq1_qcd::tr8(fenv));
      const auto _interp49 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp50 = ntRe(DiFfRG::zaaqbq1_qcd::tr18(fenv));
      const auto _interp51 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp52 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp53 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp54 = zq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp55 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp56 = zq.at(_ix4);
      const auto _interp57 = ntRe(DiFfRG::zaaqbq1_qcd::tr11(fenv));
      const auto _interp58 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp59 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp60 = zq.at(_ix5);
      const auto _interp61 = ntRe(DiFfRG::zaaqbq1_qcd::tr10(fenv));
      const auto _interp62 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp63 = ntRe(DiFfRG::zaaqbq1_qcd::tr13(fenv));
      const auto _interp64 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp65 = zq.at(_ix7);
      const auto _interp66 = ntRe(DiFfRG::zaaqbq1_qcd::tr15(fenv));
      const auto _interp67 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp68 = ntRe(DiFfRG::zaaqbq1_qcd::tr19(fenv));
      const auto _interp69 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp26 + _interp20 * powr<-1>(l1));
      const auto _den3 = powr<-2>(_interp4 * powr<2>(l1) + _interp2 * nthk0);
      const auto _den4 = powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _den6 = powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _den7 = powr<-1>(_interp5 * nthk0 + _interp6 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den8 = powr<-1>(_interp15 * nthk0 + _interp16 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den9 = powr<-1>(_interp28 * nthk0 + _interp29 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _den10 = powr<-1>(_interp34 + _interp33 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den11 = powr<-1>(_interp47 + _interp46 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den12 = powr<-1>(_interp54 + _interp53 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den13 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den14 = powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den15 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _den16 = 3. * powr<-1>(3. * _interp7 * nthk0 + _interp8 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den17 = 3. * powr<-1>(3. * _interp22 * nthk0 + _interp23 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den18 = powr<-1>(_interp36 * nthk0 + _interp37 * (powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den19 = powr<-1>(_interp41 * nthk0 + _interp42 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _den20 = powr<-1>(_interp60 + _interp59 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den21 = powr<-1>(_interp65 + _interp64 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den22 = powr<-1>(_interp56 + _interp55 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));
      // clang-format off
      return fma(0.001630434782608696, _den17 * _den2 * _interp19 * _interp24 * _interp25 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.001630434782608696, _den10 * _den2 * _den4 * _den9 * _interp27 * _interp30 * _interp31 * _interp32 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.001630434782608696, _den10 * _den18 * _den2 * _den4 * _interp31 * _interp35 * _interp38 * _interp39 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(-0.001630434782608696, _den11 * _den19 * _den2 * _den5 * _interp40 * _interp43 * _interp44 * _interp45 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(-0.001630434782608696, _den11 * _den17 * _den2 * _den5 * _interp24 * _interp44 * _interp48 * _interp49 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(-0.000815217391304348, _den13 * _den19 * _den2 * _den20 * _interp11 * _interp45 * _interp57 * _interp58 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.000815217391304348, _den13 * _den2 * _den20 * _den9 * _interp11 * _interp32 * _interp61 * _interp62 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.000815217391304348, _den16 * _den3 * _interp11 * _interp12 * _interp13 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(-0.000815217391304348, _den12 * _den15 * _den22 * _den3 * _den6 * _interp45 * _interp50 * _interp51 * _interp52 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den14 * _den21 * _den3 * _interp24 * _interp25 * _interp63 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(-0.001630434782608696, _den12 * _den14 * _den21 * _den3 * _den6 * _interp24 * _interp52 * _interp66 * _interp67 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(-0.001630434782608696, _den14 * _den15 * _den21 * _den22 * _den3 * _interp24 * _interp45 * _interp68 * _interp69 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den16 * _den3 * _den7 * _interp1 * _interp10 * _interp11 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.000815217391304348, _den3 * _den7 * _den8 * _interp14 * _interp17 * _interp18 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), 0.))))))))))))));
      // clang-format on
    }
    #elif NT_ZAAQBQ1_QCD_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2, const double& nthk0, const double& nthk1, const double& nthk2)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::zaaqbq1_qcd::nenv) > 0 ? (DiFfRG::zaaqbq1_qcd::nenv) : 1];
      DiFfRG::zaaqbq1_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp2 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp3 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix5 = ZA.index(l1);
      const auto _interp4 = ZA.at(_ix5);
      const auto _interp7 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp8 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp11 = ZAAqbq1(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _ix2 = ZAAqbq1.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp24 = ZAAqbq1.at(_ix2);
      const auto _interp25 = ZAAqbq2.at(_ix2);
      const auto _interp31 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp36 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp40 = ZAAqbq1(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p3, l1 * p, fma(2., powr<2>(p), powr<2>(l1))))));
      const auto _interp42 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp48 = ZAAqbq1(0.7071067811865475 * sqrt(fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp49 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp51 = ZAAqbq1(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(2., powr<2>(p), powr<2>(l1)))))));
      const auto _interp56 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den3 = powr<-2>(fma(_interp4, powr<2>(l1), fma(_interp2, nthk0, 0.)));
      const auto _den21 = powr<-1>(fma(_interp7, nthk0, fma(_interp8, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _interp12 = ntRe(DiFfRG::zaaqbq1_qcd::tr2(fenv));
      const auto _interp13 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp19 = ntRe(DiFfRG::zaaqbq1_qcd::tr3(fenv));
      const auto _interp20 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp21 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp22 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix8 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp23 = ZA.at(_ix8);
      const auto _interp26 = zq.at(_ix5);
      const auto _interp27 = ntIm(DiFfRG::zaaqbq1_qcd::tr4(fenv));
      const auto _interp28 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _ix11 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp29 = ZA.at(_ix11);
      const auto _ix1 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp30 = ZA3.at(_ix1);
      const auto _interp32 = ntIm(DiFfRG::zaaqbq1_qcd::tr5(fenv));
      const auto _interp33 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _ix6 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp34 = ZA.at(_ix6);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp35 = ZA3.at(_ix0);
      const auto _interp37 = ntRe(DiFfRG::zaaqbq1_qcd::tr6(fenv));
      const auto _interp38 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _ix12 = ZA.index(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp39 = ZA.at(_ix12);
      const auto _ix3 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp41 = ZAqbq1.at(_ix3);
      const auto _interp43 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix7 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp44 = zq.at(_ix7);
      const auto _interp45 = ntRe(DiFfRG::zaaqbq1_qcd::tr7(fenv));
      const auto _interp46 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix9 = ZA.index(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp47 = ZA.at(_ix9);
      const auto _interp50 = ntRe(DiFfRG::zaaqbq1_qcd::tr9(fenv));
      const auto _ix4 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp52 = ZAqbq1.at(_ix4);
      const auto _interp53 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix10 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp54 = zq.at(_ix10);
      const auto _interp55 = ntRe(DiFfRG::zaaqbq1_qcd::tr8(fenv));
      const auto _interp70 = ntRe(DiFfRG::zaaqbq1_qcd::tr11(fenv));
      const auto _interp71 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp72 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp73 = zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp74 = ntRe(DiFfRG::zaaqbq1_qcd::tr10(fenv));
      const auto _interp75 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _den13 = powr<-1>(_interp44 + _interp43 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den14 = powr<-1>(_interp54 + _interp53 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den17 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den22 = 3. * powr<-1>(3. * _interp22 * nthk0 + _interp23 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den23 = powr<-1>(_interp46 * nthk0 + _interp47 * (powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den25 = powr<-1>(_interp73 + _interp72 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp1 = ntRe(DiFfRG::zaaqbq1_qcd::tr1(fenv));
      const auto _interp5 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp6 = ZA.at(_ix7);
      const auto _interp9 = ZA3.at(_ix3);
      const auto _interp10 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp14 = ntRe(DiFfRG::zaaqbq1_qcd::tr0(fenv));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp16 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp17 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp18 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp57 = ntIm(DiFfRG::zaaqbq1_qcd::tr17(fenv));
      const auto _interp58 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp59 = ZA.at(_ix10);
      const auto _interp60 = ZA3.at(_ix4);
      const auto _interp61 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp62 = zq.at(_ix6);
      const auto _interp63 = ntRe(DiFfRG::zaaqbq1_qcd::tr18(fenv));
      const auto _interp64 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp65 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp66 = zq.at(_ix11);
      const auto _interp67 = ntIm(DiFfRG::zaaqbq1_qcd::tr12(fenv));
      const auto _interp68 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp69 = zq.at(_ix12);
      const auto _interp76 = ntRe(DiFfRG::zaaqbq1_qcd::tr13(fenv));
      const auto _interp77 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp78 = zq.at(_ix8);
      const auto _interp79 = ntIm(DiFfRG::zaaqbq1_qcd::tr14(fenv));
      const auto _interp80 = ntRe(DiFfRG::zaaqbq1_qcd::tr15(fenv));
      const auto _interp81 = ZAqbq1.at(_ix1);
      const auto _interp82 = ntRe(DiFfRG::zaaqbq1_qcd::tr19(fenv));
      const auto _interp83 = ZAqbq1.at(_ix0);
      const auto _interp84 = ntIm(DiFfRG::zaaqbq1_qcd::tr16(fenv));
      const auto _interp85 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp86 = zq.at(_ix9);
      const auto _den15 = powr<-1>(_interp66 + _interp65 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den16 = powr<-1>(_interp69 + _interp68 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den18 = powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den19 = 3. * powr<-1>(3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p));
      const auto _den26 = powr<-1>(_interp78 + _interp77 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den27 = powr<-1>(_interp86 + _interp85 * sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den28 = powr<-1>(_interp62 + _interp61 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den21 + _den22 + _den23 + _den25 + _den26 + _den27 + _den28 + _den3 + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den2 = powr<-2>(_interp26 + _interp20 * powr<-1>(l1));
        const auto _den4 = powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den11 = powr<-1>(_interp28 * nthk0 + _interp29 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
        const auto _den12 = powr<-1>(_interp38 * nthk0 + _interp39 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
        const auto _den24 = powr<-1>(_interp33 * nthk0 + _interp34 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(0.001630434782608696, _den2 * _den22 * _interp19 * _interp24 * _interp25 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.001630434782608696, _den11 * _den2 * _den22 * _interp24 * _interp27 * _interp30 * _interp31 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.001630434782608696, _den2 * _den22 * _den24 * _interp24 * _interp32 * _interp35 * _interp36 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.001630434782608696, _den12 * _den13 * _den2 * _den4 * _interp37 * _interp40 * _interp41 * _interp42 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.001630434782608696, _den13 * _den2 * _den23 * _den4 * _interp41 * _interp45 * _interp48 * _interp49 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(-0.001630434782608696, _den14 * _den2 * _den24 * _den5 * _interp36 * _interp50 * _interp51 * _interp52 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(-0.001630434782608696, _den14 * _den2 * _den22 * _den5 * _interp24 * _interp52 * _interp55 * _interp56 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(-0.000815217391304348, _den17 * _den2 * _den24 * _den25 * _interp11 * _interp36 * _interp70 * _interp71 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.000815217391304348, _den12 * _den17 * _den2 * _den25 * _interp11 * _interp42 * _interp74 * _interp75 * (_interp21 - _interp20 * etaQ) * powr<-5>(l1), fma(0.000815217391304348, _den21 * _den3 * _interp11 * _interp12 * _interp13 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), 0.))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den6 = powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den7 = powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den8 = powr<-1>(_interp5 * nthk0 + _interp6 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den9 = powr<-1>(_interp58 * nthk0 + _interp59 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp15 * nthk0 + _interp16 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den20 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
        // clang-format off
        _acc += fma(0.001630434782608696, _den20 * _den28 * _den3 * _den9 * _interp36 * _interp51 * _interp57 * _interp60 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(-0.000815217391304348, _den15 * _den20 * _den28 * _den3 * _den6 * _interp31 * _interp36 * _interp63 * _interp64 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den18 * _den26 * _den3 * _interp24 * _interp25 * _interp76 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den18 * _den26 * _den3 * _den9 * _interp24 * _interp56 * _interp60 * _interp79 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(-0.001630434782608696, _den15 * _den18 * _den26 * _den3 * _den6 * _interp24 * _interp31 * _interp80 * _interp81 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(-0.001630434782608696, _den18 * _den20 * _den26 * _den28 * _den3 * _interp24 * _interp36 * _interp82 * _interp83 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den21 * _den3 * _den8 * _interp1 * _interp10 * _interp11 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.000815217391304348, _den10 * _den3 * _den8 * _interp14 * _interp17 * _interp18 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den16 * _den3 * _den7 * _den8 * _interp40 * _interp42 * _interp67 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), fma(0.001630434782608696, _den19 * _den27 * _den3 * _den8 * _interp48 * _interp49 * _interp84 * _interp9 * (_interp3 * nthk0 + _interp2 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), 0.))))))))));
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
      double fenv[(DiFfRG::zaaqbq1_qcd::nenv) > 0 ? (DiFfRG::zaaqbq1_qcd::nenv) : 1];
      DiFfRG::zaaqbq1_qcd::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix5 = ZA.index(l1);
      const auto _interp3 = ZA.at(_ix5);
      const auto _interp6 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp7 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp10 = ZAAqbq1(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _ix2 = ZAAqbq1.index(0.7071067811865475 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp20 = ZAAqbq1.at(_ix2);
      const auto _interp21 = ZAAqbq2.at(_ix2);
      const auto _interp26 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp30 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp33 = ZAAqbq1(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p3, l1 * p, fma(2., powr<2>(p), powr<2>(l1))))));
      const auto _interp35 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp40 = ZAAqbq1(0.7071067811865475 * sqrt(fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp41 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp42 = ZAAqbq1(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(2., powr<2>(p), powr<2>(l1)))))));
      const auto _interp46 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den3 = powr<-2>(fma(_interp3, powr<2>(l1), fma(_interp1, nthk0, 0.)));
      const auto _den21 = powr<-1>(fma(_interp6, nthk0, fma(_interp7, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _interp11 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp16 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp17 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix8 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp19 = ZA.at(_ix8);
      const auto _interp22 = zq.at(_ix5);
      const auto _interp23 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _ix11 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp24 = ZA.at(_ix11);
      const auto _ix1 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp25 = ZA3.at(_ix1);
      const auto _interp27 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _ix6 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp28 = ZA.at(_ix6);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp29 = ZA3.at(_ix0);
      const auto _interp31 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _ix12 = ZA.index(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp32 = ZA.at(_ix12);
      const auto _ix3 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp34 = ZAqbq1.at(_ix3);
      const auto _interp36 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix7 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp37 = zq.at(_ix7);
      const auto _interp38 = RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix9 = ZA.index(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp39 = ZA.at(_ix9);
      const auto _ix4 = ZAqbq1.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp43 = ZAqbq1.at(_ix4);
      const auto _interp44 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _ix10 = zq.index(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp45 = zq.at(_ix10);
      const auto _interp57 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp58 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp59 = zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp60 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _den13 = powr<-1>(_interp37 + _interp36 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))));
      const auto _den14 = powr<-1>(_interp45 + _interp44 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den17 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den22 = 3. * powr<-1>(3. * _interp18 * nthk0 + _interp19 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den23 = powr<-1>(_interp38 * nthk0 + _interp39 * (powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den25 = powr<-1>(_interp59 + _interp58 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _interp5 = ZA.at(_ix7);
      const auto _interp8 = ZA3.at(_ix3);
      const auto _interp9 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp12 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp13 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp14 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp15 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp47 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp48 = ZA.at(_ix10);
      const auto _interp49 = ZA3.at(_ix4);
      const auto _interp50 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _interp51 = zq.at(_ix6);
      const auto _interp52 = ZAAqbq1(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp53 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp54 = zq.at(_ix11);
      const auto _interp55 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp56 = zq.at(_ix12);
      const auto _interp61 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp62 = zq.at(_ix8);
      const auto _interp63 = ZAqbq1.at(_ix1);
      const auto _interp64 = ZAqbq1.at(_ix0);
      const auto _interp65 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp66 = zq.at(_ix9);
      const auto _den15 = powr<-1>(_interp54 + _interp53 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den16 = powr<-1>(_interp56 + _interp55 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den18 = powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _den19 = 3. * powr<-1>(3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p));
      const auto _den26 = powr<-1>(_interp62 + _interp61 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den27 = powr<-1>(_interp66 + _interp65 * sqrt(powr<-1>(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den28 = powr<-1>(_interp51 + _interp50 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p))));
      // clang-format off
      using _T = decltype(complex<double>(0.,-0.001630434782608696) + _den1 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den21 + _den22 + _den23 + _den25 + _den26 + _den27 + _den28 + _den3 + DiFfRG::zaaqbq1_qcd::tr0(fenv) + DiFfRG::zaaqbq1_qcd::tr10(fenv) + DiFfRG::zaaqbq1_qcd::tr11(fenv) + DiFfRG::zaaqbq1_qcd::tr12(fenv) + DiFfRG::zaaqbq1_qcd::tr13(fenv) + DiFfRG::zaaqbq1_qcd::tr14(fenv) + DiFfRG::zaaqbq1_qcd::tr15(fenv) + DiFfRG::zaaqbq1_qcd::tr16(fenv) + DiFfRG::zaaqbq1_qcd::tr17(fenv) + DiFfRG::zaaqbq1_qcd::tr18(fenv) + DiFfRG::zaaqbq1_qcd::tr19(fenv) + DiFfRG::zaaqbq1_qcd::tr1(fenv) + DiFfRG::zaaqbq1_qcd::tr2(fenv) + DiFfRG::zaaqbq1_qcd::tr3(fenv) + DiFfRG::zaaqbq1_qcd::tr4(fenv) + DiFfRG::zaaqbq1_qcd::tr5(fenv) + DiFfRG::zaaqbq1_qcd::tr6(fenv) + DiFfRG::zaaqbq1_qcd::tr7(fenv) + DiFfRG::zaaqbq1_qcd::tr8(fenv) + DiFfRG::zaaqbq1_qcd::tr9(fenv) + _interp1 + _interp10 + _interp11 + _interp12 + _interp13 + _interp14 + _interp15 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp7 + _interp8 + _interp9 + cosl1p1 + cosl1p2 + cosl1p3 + etaQ + k + l1 + nthk0 + nthk1 + nthk2 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den2 = powr<-2>(_interp22 + _interp16 * powr<-1>(l1));
        const auto _den4 = powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den11 = powr<-1>(_interp23 * nthk0 + _interp24 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
        const auto _den12 = powr<-1>(_interp31 * nthk0 + _interp32 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
        const auto _den24 = powr<-1>(_interp27 * nthk0 + _interp28 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
        const auto _cse1_k1 = powr<-5>(l1);
        const auto _cse2_k1 = -_interp16 * etaQ;
        const auto _cse3_k1 = _cse2_k1 + _interp17;
        // clang-format off
        _acc += fma(0.001630434782608696, _cse1_k1 * _cse3_k1 * _den2 * _den22 * DiFfRG::zaaqbq1_qcd::tr3(fenv) * _interp20 * _interp21, fma(complex<double>(0.,-0.001630434782608696), _cse1_k1 * _cse3_k1 * _den11 * _den2 * _den22 * DiFfRG::zaaqbq1_qcd::tr4(fenv) * _interp20 * _interp25 * _interp26, fma(complex<double>(0.,-0.001630434782608696), _cse1_k1 * _cse3_k1 * _den2 * _den22 * _den24 * DiFfRG::zaaqbq1_qcd::tr5(fenv) * _interp20 * _interp29 * _interp30, fma(0.001630434782608696, _cse1_k1 * _cse3_k1 * _den12 * _den13 * _den2 * _den4 * DiFfRG::zaaqbq1_qcd::tr6(fenv) * _interp33 * _interp34 * _interp35, fma(0.001630434782608696, _cse1_k1 * _cse3_k1 * _den13 * _den2 * _den23 * _den4 * DiFfRG::zaaqbq1_qcd::tr7(fenv) * _interp34 * _interp40 * _interp41, fma(-0.001630434782608696, _cse1_k1 * _cse3_k1 * _den14 * _den2 * _den24 * _den5 * DiFfRG::zaaqbq1_qcd::tr9(fenv) * _interp30 * _interp42 * _interp43, fma(-0.001630434782608696, _cse1_k1 * _cse3_k1 * _den14 * _den2 * _den22 * _den5 * DiFfRG::zaaqbq1_qcd::tr8(fenv) * _interp20 * _interp43 * _interp46, fma(-0.000815217391304348, _cse1_k1 * _cse3_k1 * _den17 * _den2 * _den24 * _den25 * DiFfRG::zaaqbq1_qcd::tr11(fenv) * _interp10 * _interp30 * _interp57, fma(0.000815217391304348, _cse1_k1 * _cse3_k1 * _den12 * _den17 * _den2 * _den25 * DiFfRG::zaaqbq1_qcd::tr10(fenv) * _interp10 * _interp35 * _interp60, fma(0.000815217391304348, _den21 * _den3 * DiFfRG::zaaqbq1_qcd::tr2(fenv) * _interp10 * _interp11 * (_interp2 * nthk0 + _interp1 * (nthk1 + 50. * _den1 * powr<6>(k) * (-nthk0 + nthk2))), 0.))))))))));
        // clang-format on
      }
      { // subkernel 2
        const auto _den6 = powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den7 = powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p));
        const auto _den8 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den9 = powr<-1>(_interp47 * nthk0 + _interp48 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp12 * nthk0 + _interp13 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den20 = powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
        const auto _cse1_k2 = _interp2 * nthk0;
        const auto _cse2_k2 = -nthk0;
        const auto _cse3_k2 = _cse2_k2 + nthk2;
        const auto _cse4_k2 = 50. * _cse3_k2 * _den1 * powr<6>(k);
        const auto _cse5_k2 = _cse4_k2 + nthk1;
        const auto _cse6_k2 = _cse5_k2 * _interp1;
        const auto _cse7_k2 = _cse1_k2 + _cse6_k2;
        // clang-format off
        _acc += 0.000815217391304348 * fma(2., _cse7_k2 * _den18 * _den26 * _den3 * DiFfRG::zaaqbq1_qcd::tr13(fenv) * _interp20 * _interp21, fma(complex<double>(0.,-2.), _cse7_k2 * _den20 * _den28 * _den3 * _den9 * DiFfRG::zaaqbq1_qcd::tr17(fenv) * _interp30 * _interp42 * _interp49, fma(complex<double>(0.,-2.), _cse7_k2 * _den18 * _den26 * _den3 * _den9 * DiFfRG::zaaqbq1_qcd::tr14(fenv) * _interp20 * _interp46 * _interp49, fma(-1., _cse7_k2 * _den15 * _den20 * _den28 * _den3 * _den6 * DiFfRG::zaaqbq1_qcd::tr18(fenv) * _interp26 * _interp30 * _interp52, fma(-2., _cse7_k2 * _den15 * _den18 * _den26 * _den3 * _den6 * DiFfRG::zaaqbq1_qcd::tr15(fenv) * _interp20 * _interp26 * _interp63, fma(-2., _cse7_k2 * _den18 * _den20 * _den26 * _den28 * _den3 * DiFfRG::zaaqbq1_qcd::tr19(fenv) * _interp20 * _interp30 * _interp64, fma(_cse7_k2, _den10 * _den3 * _den8 * DiFfRG::zaaqbq1_qcd::tr0(fenv) * _interp14 * _interp15 * _interp8, fma(complex<double>(0.,-2.), _cse7_k2 * _den16 * _den3 * _den7 * _den8 * DiFfRG::zaaqbq1_qcd::tr12(fenv) * _interp33 * _interp35 * _interp8, fma(complex<double>(0.,-2.), _cse7_k2 * _den19 * _den27 * _den3 * _den8 * DiFfRG::zaaqbq1_qcd::tr16(fenv) * _interp40 * _interp41 * _interp8, fma(2., _cse7_k2 * _den21 * _den3 * _den8 * DiFfRG::zaaqbq1_qcd::tr1(fenv) * _interp10 * _interp8 * _interp9, 0.))))))))));
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
      const auto _interp1 = ZAAqbq1(p);
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
using DiFfRG::ZAAqbq1_kernel;