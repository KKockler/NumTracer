#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "ZA4_gpu_kernels.hh"
#include "numtrace_verdict.hh"

namespace DiFfRG
{
  class ZA4_gpu_kernel
  {
    public:
    #if NT_ZA4_GPU_VERDICT == 2   // Pure: the Complex -> Re projection is exact
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::za4_gpu::nenv) > 0 ? (DiFfRG::za4_gpu::nenv) : 1];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_4 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_6 = Mq(l1);
      const double dr_7 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      const double dr_8 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_10 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_13 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_14 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_15 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_16 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_17 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_20 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_21 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_22 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_23 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      DiFfRG::za4_gpu::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23);
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp9 = ZA(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp10 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp11 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp12 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp13 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp14 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp15 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp16 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp17 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp19 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp20 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp21 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp22 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp24 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp25 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp26 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp29 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp30 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp31 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp34 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp35 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp36 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp37 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp39 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp41 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp42 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp43 = ZA3(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp46 = RB(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp47 = ZA(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp48 = ZA3(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp52 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp53 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp54 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp57 = RB(powr<2>(k), fma(2., cosl1p1 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp58 = ZA(sqrt(fma(2., cosl1p1 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp59 = ZA3(0.816496580927726 * sqrt(fma(cosl1p1, l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp96 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp97 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp100 = Zc(k);
      const auto _interp101 = dtZc(k);
      const auto _interp102 = Zc(1.02 * k);
      const auto _interp103 = Zc(l1);
      const auto _interp104 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp105 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp106 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp109 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp110 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp111 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp114 = ZAcbc(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp115 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp116 = Zc(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp122 = ZAcbc(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp125 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp128 = ZAcbc(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp131 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp133 = Zc(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp136 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp138 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp141 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp3, _interp4, fma(_interp7, powr<2>(l1), 0.)));
      const auto _den3 = powr<-2>(fma(-1., _interp100 * _interp4, fma(-1., _interp103 * powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp24, _interp3, fma(_interp25, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp29, _interp3, fma(_interp30, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den6 = powr<-1>(fma(_interp3, _interp41, fma(_interp42, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den7 = powr<-1>(fma(_interp3, _interp8, fma(_interp9, powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(_interp3, _interp46, fma(_interp47, powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den9 = powr<-1>(fma(-1., _interp100 * _interp24, fma(-1., _interp104 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)), 0.)));
      const auto _den10 = powr<-1>(fma(-1., _interp100 * _interp29, fma(-1., _interp111 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)), 0.)));
      const auto _den12 = powr<-1>(fma(-1., _interp100 * _interp8, fma(-1., _interp133 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)), 0.)));
      const auto _den13 = powr<-1>(fma(-1., _interp100 * _interp46, fma(-1., _interp116 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)), 0.)));
      const auto _den14 = powr<-1>(fma(_interp3, _interp34, fma(_interp35, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den15 = powr<-1>(fma(_interp19, _interp3, fma(_interp20, powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den16 = powr<-1>(fma(_interp12, _interp3, fma(_interp13, powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den17 = powr<-1>(fma(_interp3, _interp57, fma(_interp58, powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den18 = powr<-1>(fma(_interp3, _interp52, fma(_interp53, powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den19 = powr<-1>(fma(_interp10, _interp3, fma(_interp11, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
      const auto _den20 = powr<-1>(fma(-1., _interp100 * _interp34, fma(-1., _interp106 * powr<2>(l1) - 2. * _interp106 * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * _interp106 * powr<2>(p), 0.)));
      const auto _den21 = powr<-1>(fma(-1., _interp100 * _interp19, fma(-1., _interp125 * powr<2>(l1) - 2. * _interp125 * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * _interp125 * powr<2>(p), 0.)));
      const auto _den22 = powr<-1>(fma(-1., _interp100 * _interp12, fma(-1., _interp138 * powr<2>(l1) - 2. * _interp138 * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * _interp138 * powr<2>(p), 0.)));
      const auto _den25 = powr<-1>(fma(-1., _interp10 * _interp100, fma(-1., _interp105 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)), 0.)));
      const auto _interp118 = ntRe(DiFfRG::za4_gpu::tr35(fenv));
      const auto _interp119 = ntRe(DiFfRG::za4_gpu::tr36(fenv));
      const auto _interp120 = ntRe(DiFfRG::za4_gpu::tr37(fenv));
      const auto _interp121 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp123 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp124 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp126 = ntRe(DiFfRG::za4_gpu::tr38(fenv));
      const auto _interp130 = ntRe(DiFfRG::za4_gpu::tr51(fenv));
      const auto _interp132 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp134 = ntRe(DiFfRG::za4_gpu::tr52(fenv));
      const auto _interp135 = ntRe(DiFfRG::za4_gpu::tr53(fenv));
      const auto _interp137 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp139 = ntRe(DiFfRG::za4_gpu::tr54(fenv));
      const auto _interp140 = ntRe(DiFfRG::za4_gpu::tr55(fenv));
      const auto _interp143 = ntRe(DiFfRG::za4_gpu::tr39(fenv));
      const auto _interp144 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp145 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp146 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp147 = ntRe(DiFfRG::za4_gpu::tr40(fenv));
      const auto _interp148 = ntRe(DiFfRG::za4_gpu::tr41(fenv));
      const auto _interp149 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp150 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp151 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp152 = ntRe(DiFfRG::za4_gpu::tr42(fenv));
      const auto _interp153 = ntRe(DiFfRG::za4_gpu::tr43(fenv));
      const auto _interp154 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp155 = ntRe(DiFfRG::za4_gpu::tr44(fenv));
      const auto _den23 = -3. * powr<-1>(3. * _interp100 * _interp57 + _interp146 * (3. * powr<2>(l1) + 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den24 = -3. * powr<-1>(3. * _interp100 * _interp52 + _interp151 * (3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _interp1 = ntRe(DiFfRG::za4_gpu::tr10(fenv));
      const auto _interp18 = ntRe(DiFfRG::za4_gpu::tr8(fenv));
      const auto _interp23 = ntRe(DiFfRG::za4_gpu::tr9(fenv));
      const auto _interp28 = ntRe(DiFfRG::za4_gpu::tr11(fenv));
      const auto _interp33 = ntRe(DiFfRG::za4_gpu::tr6(fenv));
      const auto _interp38 = ntRe(DiFfRG::za4_gpu::tr7(fenv));
      const auto _interp95 = ntRe(DiFfRG::za4_gpu::tr45(fenv));
      const auto _interp98 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp99 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp107 = ntRe(DiFfRG::za4_gpu::tr46(fenv));
      const auto _interp108 = ntRe(DiFfRG::za4_gpu::tr47(fenv));
      const auto _interp112 = ntRe(DiFfRG::za4_gpu::tr48(fenv));
      const auto _interp113 = ntRe(DiFfRG::za4_gpu::tr33(fenv));
      const auto _interp117 = ntRe(DiFfRG::za4_gpu::tr34(fenv));
      const auto _interp127 = ntRe(DiFfRG::za4_gpu::tr49(fenv));
      const auto _interp129 = ntRe(DiFfRG::za4_gpu::tr50(fenv));
      const auto _interp142 = ntRe(DiFfRG::za4_gpu::tr56(fenv));
      const auto _interp40 = ntRe(DiFfRG::za4_gpu::tr2(fenv));
      const auto _interp44 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp45 = ntRe(DiFfRG::za4_gpu::tr0(fenv));
      const auto _interp49 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp50 = ntRe(DiFfRG::za4_gpu::tr1(fenv));
      const auto _interp51 = ntRe(DiFfRG::za4_gpu::tr4(fenv));
      const auto _interp55 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp56 = ntRe(DiFfRG::za4_gpu::tr3(fenv));
      const auto _interp60 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp61 = ntRe(DiFfRG::za4_gpu::tr5(fenv));
      const auto _interp62 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp63 = ntRe(DiFfRG::za4_gpu::tr24(fenv));
      const auto _interp64 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp65 = ntRe(DiFfRG::za4_gpu::tr25(fenv));
      const auto _interp66 = ntRe(DiFfRG::za4_gpu::tr13(fenv));
      const auto _interp67 = ntRe(DiFfRG::za4_gpu::tr12(fenv));
      const auto _interp68 = ntRe(DiFfRG::za4_gpu::tr26(fenv));
      const auto _interp69 = ntRe(DiFfRG::za4_gpu::tr27(fenv));
      const auto _interp70 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp71 = ntRe(DiFfRG::za4_gpu::tr14(fenv));
      const auto _interp72 = ntRe(DiFfRG::za4_gpu::tr28(fenv));
      const auto _interp73 = ntRe(DiFfRG::za4_gpu::tr29(fenv));
      const auto _interp74 = ntRe(DiFfRG::za4_gpu::tr16(fenv));
      const auto _interp75 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp76 = ntRe(DiFfRG::za4_gpu::tr30(fenv));
      const auto _interp77 = ntRe(DiFfRG::za4_gpu::tr31(fenv));
      const auto _interp78 = ntRe(DiFfRG::za4_gpu::tr32(fenv));
      const auto _interp79 = ntRe(DiFfRG::za4_gpu::tr15(fenv));
      const auto _interp80 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p1 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp81 = ntRe(DiFfRG::za4_gpu::tr17(fenv));
      const auto _interp82 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp83 = ntRe(DiFfRG::za4_gpu::tr18(fenv));
      const auto _interp84 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp85 = ntRe(DiFfRG::za4_gpu::tr19(fenv));
      const auto _interp86 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp87 = ntRe(DiFfRG::za4_gpu::tr20(fenv));
      const auto _interp88 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp89 = ntRe(DiFfRG::za4_gpu::tr21(fenv));
      const auto _interp90 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (2. * cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp91 = ntRe(DiFfRG::za4_gpu::tr22(fenv));
      const auto _interp92 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp93 = ntRe(DiFfRG::za4_gpu::tr23(fenv));
      const auto _interp94 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den12 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den2 + _den20 + _den21 + _den22 + _den23 + _den24 + _den25 + _den3 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp122 + _interp123 + _interp124 + _interp125 + _interp126 + _interp127 + _interp128 + _interp129 + _interp13 + _interp130 + _interp131 + _interp132 + _interp133 + _interp134 + _interp135 + _interp136 + _interp137 + _interp138 + _interp139 + _interp14 + _interp140 + _interp141 + _interp142 + _interp143 + _interp144 + _interp145 + _interp146 + _interp147 + _interp148 + _interp149 + _interp15 + _interp150 + _interp151 + _interp152 + _interp153 + _interp154 + _interp155 + _interp16 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p2 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den11 = -powr<-1>(_interp100 * _interp41 + _interp124 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(-0.0004251700680272108, _den10 * _den13 * _den20 * _den3 * _interp109 * _interp110 * _interp114 * _interp115 * _interp118 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(0.0004251700680272108, _den10 * _den13 * _den20 * _den3 * _interp109 * _interp110 * _interp114 * _interp115 * _interp119 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(-0.0004251700680272108, _den10 * _den13 * _den23 * _den3 * _interp109 * _interp114 * _interp143 * _interp144 * _interp145 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(0.0004251700680272108, _den10 * _den13 * _den23 * _den3 * _interp109 * _interp114 * _interp144 * _interp145 * _interp147 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(0.0004251700680272108, _den11 * _den21 * _den3 * _den9 * _interp120 * _interp121 * _interp122 * _interp123 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(-0.0004251700680272108, _den11 * _den21 * _den3 * _den9 * _interp121 * _interp122 * _interp123 * _interp126 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(-0.0004251700680272108, _den11 * _den24 * _den3 * _den9 * _interp121 * _interp148 * _interp149 * _interp150 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(0.0004251700680272108, _den11 * _den24 * _den3 * _den9 * _interp121 * _interp149 * _interp150 * _interp152 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(-0.0004251700680272108, _den13 * _den24 * _den3 * _den9 * _interp114 * _interp149 * _interp153 * _interp154 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(0.0004251700680272108, _den13 * _den24 * _den3 * _den9 * _interp114 * _interp149 * _interp154 * _interp155 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(0.0004251700680272108, _den12 * _den21 * _den25 * _den3 * _interp128 * _interp130 * _interp131 * _interp132 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(-0.0004251700680272108, _den12 * _den21 * _den25 * _den3 * _interp128 * _interp131 * _interp132 * _interp134 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(0.0004251700680272108, _den10 * _den22 * _den25 * _den3 * _interp109 * _interp135 * _interp136 * _interp137 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(-0.0004251700680272108, _den10 * _den22 * _den25 * _den3 * _interp109 * _interp136 * _interp137 * _interp139 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(0.0004251700680272108, _den12 * _den22 * _den25 * _den3 * _interp131 * _interp136 * _interp140 * _interp141 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, 0.)))))))))))))));
        // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(-0.0004251700680272108, _den12 * _den22 * _den25 * _den3 * _interp131 * _interp136 * _interp141 * _interp142 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(0.0004251700680272108, _den21 * _den25 * _den3 * _den9 * _interp122 * _interp127 * _interp128 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp97, fma(-0.0004251700680272108, _den21 * _den25 * _den3 * _den9 * _interp122 * _interp128 * _interp129 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp97, fma(-0.0004251700680272108, _den13 * _den20 * _den3 * _den9 * _interp113 * _interp114 * _interp115 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp98, fma(0.0004251700680272108, _den13 * _den20 * _den3 * _den9 * _interp114 * _interp115 * _interp117 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp98, fma(0.0004251700680272108, _den10 * _den20 * _den25 * _den3 * _interp108 * _interp109 * _interp110 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97 * _interp99, fma(-0.0004251700680272108, _den10 * _den20 * _den25 * _den3 * _interp109 * _interp110 * _interp112 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97 * _interp99, fma(-0.0004251700680272108, _den20 * _den25 * _den3 * _den9 * _interp107 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp97 * _interp98 * _interp99, fma(0.0004251700680272108, _den20 * _den25 * _den3 * _den9 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp95 * _interp96 * _interp97 * _interp98 * _interp99, fma(0.0004251700680272108, _den16 * _den19 * _den2 * _den7 * _interp1 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den19 * _den2 * _den7 * _interp14 * _interp15 * _interp18 * _interp21 * _interp22 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den19 * _den2 * _den4 * _interp15 * _interp22 * _interp23 * _interp26 * _interp27 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den19 * _den2 * _den5 * _interp15 * _interp16 * _interp28 * _interp31 * _interp32 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den19 * _den2 * _den5 * _interp15 * _interp31 * _interp33 * _interp36 * _interp37 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den19 * _den2 * _den4 * _interp15 * _interp26 * _interp37 * _interp38 * _interp39 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.)))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.0004251700680272108, _den15 * _den2 * _den4 * _den6 * _interp26 * _interp27 * _interp40 * _interp43 * _interp44 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den5 * _den8 * _interp31 * _interp36 * _interp45 * _interp48 * _interp49 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den4 * _den8 * _interp26 * _interp39 * _interp48 * _interp49 * _interp50 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den4 * _den6 * _interp26 * _interp43 * _interp51 * _interp54 * _interp55 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den5 * _den8 * _interp31 * _interp48 * _interp56 * _interp59 * _interp60 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den4 * _den8 * _interp26 * _interp48 * _interp54 * _interp61 * _interp62 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den5 * _interp31 * _interp36 * _interp63 * _interp64 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den4 * _interp26 * _interp39 * _interp64 * _interp65 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den19 * _den2 * _interp15 * _interp37 * _interp64 * _interp66 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * _den8 * _interp48 * _interp49 * _interp64 * _interp67 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den14 * _den2 * powr<2>(_interp64) * _interp68 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den7 * _interp14 * _interp21 * _interp69 * _interp70 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den19 * _den2 * _interp15 * _interp22 * _interp70 * _interp71 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den4 * _interp26 * _interp27 * _interp70 * _interp72 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * powr<2>(_interp70) * _interp73 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.)))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(0.0004251700680272108, _den16 * _den19 * _den2 * _interp15 * _interp16 * _interp74 * _interp75 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den7 * _interp14 * _interp17 * _interp75 * _interp76 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den5 * _interp31 * _interp32 * _interp75 * _interp77 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * powr<2>(_interp75) * _interp78 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den5 * _interp31 * _interp59 * _interp79 * _interp80 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den4 * _interp26 * _interp54 * _interp81 * _interp82 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den6 * _interp26 * _interp43 * _interp83 * _interp84 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den7 * _interp14 * _interp15 * _interp85 * _interp86 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den5 * _interp15 * _interp31 * _interp87 * _interp88 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * _interp15 * _interp26 * _interp89 * _interp90 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den8 * _interp26 * _interp48 * _interp91 * _interp92 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den5 * _den8 * _interp31 * _interp48 * _interp93 * _interp94 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.))))))))))));
      // clang-format on
      }
      return _acc;
    }
    #elif NT_ZA4_GPU_VERDICT == 1   // RePart: real value via complex trace(s), re/im split
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::za4_gpu::nenv) > 0 ? (DiFfRG::za4_gpu::nenv) : 1];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_4 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_6 = Mq(l1);
      const double dr_7 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      const double dr_8 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_10 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_13 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_14 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_15 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_16 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_17 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_20 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_21 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_22 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_23 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      DiFfRG::za4_gpu::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23);
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp3 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp5 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp7 = ZA(l1);
      const auto _interp8 = RB(powr<2>(k), fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp9 = ZA(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp10 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp11 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp12 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp13 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp14 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp15 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp16 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp17 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp19 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp20 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp21 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp22 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp24 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp25 = ZA(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp26 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp27 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp29 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp30 = ZA(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp31 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp32 = ZA3(0.816496580927726 * sqrt(fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp34 = RB(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp41 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp42 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp43 = ZA3(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp46 = RB(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp47 = ZA(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp48 = ZA3(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp52 = RB(powr<2>(k), fma(2., cosl1p2 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp53 = ZA(sqrt(fma(2., cosl1p2 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp54 = ZA3(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(cosl1p2, l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp57 = RB(powr<2>(k), fma(2., cosl1p1 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp58 = ZA(sqrt(fma(2., cosl1p1 * l1 * p, fma(2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp59 = ZA3(0.816496580927726 * sqrt(fma(cosl1p1, l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp70 = ZA4(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp96 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp97 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp100 = Zc(k);
      const auto _interp101 = dtZc(k);
      const auto _interp102 = Zc(1.02 * k);
      const auto _interp103 = Zc(l1);
      const auto _interp104 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp105 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp106 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp109 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp110 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp111 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp114 = ZAcbc(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp115 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp116 = Zc(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp122 = ZAcbc(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp125 = Zc(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp128 = ZAcbc(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp131 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp133 = Zc(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp136 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp138 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp157 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp158 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp161 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp162 = Zq(k);
      const auto _interp163 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp164 = dtZq(k);
      const auto _interp165 = Zq(1.02 * k);
      const auto _interp166 = Mq(l1);
      const auto _interp167 = Zq(l1);
      const auto _interp168 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp169 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp170 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp171 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp172 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp173 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp179 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp181 = Mq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp182 = RF(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp183 = Zq(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp186 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp188 = Mq(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp189 = RF(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp190 = Zq(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp195 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp198 = Mq(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp199 = RF(powr<2>(k), fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp200 = Zq(sqrt(fma(2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp201 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp202 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp203 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp206 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp209 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp210 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp211 = Mq(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp212 = RF(powr<2>(k), fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp213 = Zq(sqrt(fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(fma(_interp3, _interp4, fma(_interp7, powr<2>(l1), 0.)));
      const auto _den3 = powr<-2>(fma(-1., _interp100 * _interp4, fma(-1., _interp103 * powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp24, _interp3, fma(_interp25, powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), 0.)));
      const auto _den5 = powr<-1>(fma(_interp29, _interp3, fma(_interp30, powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den6 = powr<-1>(fma(_interp3, _interp41, fma(_interp42, powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), 0.)));
      const auto _den7 = powr<-1>(fma(_interp3, _interp8, fma(_interp9, powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den8 = powr<-1>(fma(_interp3, _interp46, fma(_interp47, powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), 0.)));
      const auto _den9 = powr<-1>(fma(-1., _interp100 * _interp24, fma(-1., _interp104 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)), 0.)));
      const auto _den10 = powr<-1>(fma(-1., _interp100 * _interp29, fma(-1., _interp111 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)), 0.)));
      const auto _den12 = powr<-1>(fma(-1., _interp100 * _interp8, fma(-1., _interp133 * (powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)), 0.)));
      const auto _den13 = powr<-1>(fma(-1., _interp100 * _interp46, fma(-1., _interp116 * (powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)), 0.)));
      const auto _den14 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp162 * _interp163 + _interp167 * l1), powr<2>(l1), powr<2>(_interp166)));
      const auto _den16 = powr<-1>(fma(_interp19, _interp3, fma(_interp20, powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den17 = powr<-1>(fma(_interp12, _interp3, fma(_interp13, powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den18 = powr<-1>(fma(_interp3, _interp57, fma(_interp58, powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den19 = powr<-1>(fma(_interp3, _interp52, fma(_interp53, powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den20 = powr<-1>(fma(_interp10, _interp3, fma(_interp11, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
      const auto _den21 = powr<-1>(fma(-1., _interp100 * _interp34, fma(-1., _interp106 * powr<2>(l1) - 2. * _interp106 * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * _interp106 * powr<2>(p), 0.)));
      const auto _den22 = powr<-1>(fma(-1., _interp100 * _interp19, fma(-1., _interp125 * powr<2>(l1) - 2. * _interp125 * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * _interp125 * powr<2>(p), 0.)));
      const auto _den23 = powr<-1>(fma(-1., _interp100 * _interp12, fma(-1., _interp138 * powr<2>(l1) - 2. * _interp138 * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * _interp138 * powr<2>(p), 0.)));
      const auto _den26 = powr<-1>(fma(-1., _interp10 * _interp100, fma(-1., _interp105 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)), 0.)));
      const auto _den27 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp170 + _interp162 * _interp169 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp168)));
      const auto _den28 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp183 + _interp162 * _interp182 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp181)));
      const auto _den29 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p), powr<2>(_interp200 + _interp162 * _interp199 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))), powr<2>(_interp198)));
      const auto _den30 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p), powr<2>(_interp213 + _interp162 * _interp212 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)))), powr<2>(_interp211)));
      const auto _den31 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), powr<2>(_interp190 + _interp162 * _interp189 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)))), powr<2>(_interp188)));
      const auto _den33 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp203 + _interp162 * _interp202 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp201)));
      const auto _den37 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp173 + _interp162 * _interp172 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp171)));
      const auto _interp118 = ntRe(DiFfRG::za4_gpu::tr35(fenv));
      const auto _interp119 = ntRe(DiFfRG::za4_gpu::tr36(fenv));
      const auto _interp120 = ntRe(DiFfRG::za4_gpu::tr37(fenv));
      const auto _interp121 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp123 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp124 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp126 = ntRe(DiFfRG::za4_gpu::tr38(fenv));
      const auto _interp130 = ntRe(DiFfRG::za4_gpu::tr51(fenv));
      const auto _interp132 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp134 = ntRe(DiFfRG::za4_gpu::tr52(fenv));
      const auto _interp135 = ntRe(DiFfRG::za4_gpu::tr53(fenv));
      const auto _interp137 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp139 = ntRe(DiFfRG::za4_gpu::tr54(fenv));
      const auto _interp143 = ntRe(DiFfRG::za4_gpu::tr39(fenv));
      const auto _interp144 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp145 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp146 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp147 = ntRe(DiFfRG::za4_gpu::tr40(fenv));
      const auto _interp148 = ntRe(DiFfRG::za4_gpu::tr41(fenv));
      const auto _interp149 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp150 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp151 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp152 = ntRe(DiFfRG::za4_gpu::tr42(fenv));
      const auto _interp153 = ntRe(DiFfRG::za4_gpu::tr43(fenv));
      const auto _interp154 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp155 = ntRe(DiFfRG::za4_gpu::tr44(fenv));
      const auto _den24 = -3. * powr<-1>(3. * _interp100 * _interp57 + _interp146 * (3. * powr<2>(l1) + 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den25 = -3. * powr<-1>(3. * _interp100 * _interp52 + _interp151 * (3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _interp1 = ntRe(DiFfRG::za4_gpu::tr10(fenv));
      const auto _interp18 = ntRe(DiFfRG::za4_gpu::tr8(fenv));
      const auto _interp23 = ntRe(DiFfRG::za4_gpu::tr9(fenv));
      const auto _interp28 = ntRe(DiFfRG::za4_gpu::tr11(fenv));
      const auto _interp95 = ntRe(DiFfRG::za4_gpu::tr45(fenv));
      const auto _interp98 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp99 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp107 = ntRe(DiFfRG::za4_gpu::tr46(fenv));
      const auto _interp108 = ntRe(DiFfRG::za4_gpu::tr47(fenv));
      const auto _interp112 = ntRe(DiFfRG::za4_gpu::tr48(fenv));
      const auto _interp113 = ntRe(DiFfRG::za4_gpu::tr33(fenv));
      const auto _interp117 = ntRe(DiFfRG::za4_gpu::tr34(fenv));
      const auto _interp127 = ntRe(DiFfRG::za4_gpu::tr49(fenv));
      const auto _interp129 = ntRe(DiFfRG::za4_gpu::tr50(fenv));
      const auto _interp140 = ntRe(DiFfRG::za4_gpu::tr55(fenv));
      const auto _interp141 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp142 = ntRe(DiFfRG::za4_gpu::tr56(fenv));
      const auto _interp33 = ntRe(DiFfRG::za4_gpu::tr6(fenv));
      const auto _interp35 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp36 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp37 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp38 = ntRe(DiFfRG::za4_gpu::tr7(fenv));
      const auto _interp39 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp40 = ntRe(DiFfRG::za4_gpu::tr2(fenv));
      const auto _interp44 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp45 = ntRe(DiFfRG::za4_gpu::tr0(fenv));
      const auto _interp49 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp50 = ntRe(DiFfRG::za4_gpu::tr1(fenv));
      const auto _interp51 = ntRe(DiFfRG::za4_gpu::tr4(fenv));
      const auto _interp55 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp56 = ntRe(DiFfRG::za4_gpu::tr3(fenv));
      const auto _interp60 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp61 = ntRe(DiFfRG::za4_gpu::tr5(fenv));
      const auto _interp62 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp63 = ntRe(DiFfRG::za4_gpu::tr24(fenv));
      const auto _interp64 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp65 = ntRe(DiFfRG::za4_gpu::tr25(fenv));
      const auto _interp66 = ntRe(DiFfRG::za4_gpu::tr13(fenv));
      const auto _interp67 = ntRe(DiFfRG::za4_gpu::tr12(fenv));
      const auto _interp68 = ntRe(DiFfRG::za4_gpu::tr26(fenv));
      const auto _interp69 = ntRe(DiFfRG::za4_gpu::tr27(fenv));
      const auto _den15 = 3. * powr<-1>(3. * _interp3 * _interp34 + _interp35 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _interp71 = ntRe(DiFfRG::za4_gpu::tr14(fenv));
      const auto _interp72 = ntRe(DiFfRG::za4_gpu::tr28(fenv));
      const auto _interp73 = ntRe(DiFfRG::za4_gpu::tr29(fenv));
      const auto _interp74 = ntRe(DiFfRG::za4_gpu::tr16(fenv));
      const auto _interp75 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp76 = ntRe(DiFfRG::za4_gpu::tr30(fenv));
      const auto _interp77 = ntRe(DiFfRG::za4_gpu::tr31(fenv));
      const auto _interp78 = ntRe(DiFfRG::za4_gpu::tr32(fenv));
      const auto _interp79 = ntRe(DiFfRG::za4_gpu::tr15(fenv));
      const auto _interp80 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p1 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp81 = ntRe(DiFfRG::za4_gpu::tr17(fenv));
      const auto _interp82 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp83 = ntRe(DiFfRG::za4_gpu::tr18(fenv));
      const auto _interp84 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp85 = ntRe(DiFfRG::za4_gpu::tr19(fenv));
      const auto _interp86 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp87 = ntRe(DiFfRG::za4_gpu::tr20(fenv));
      const auto _interp88 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp89 = ntRe(DiFfRG::za4_gpu::tr21(fenv));
      const auto _interp90 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (2. * cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp91 = ntRe(DiFfRG::za4_gpu::tr22(fenv));
      const auto _interp92 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp93 = ntRe(DiFfRG::za4_gpu::tr23(fenv));
      const auto _interp94 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp156 = ntIm(DiFfRG::za4_gpu::tr69(fenv));
      const auto _interp159 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp160 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp174 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp175 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp176 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp177 = ntIm(DiFfRG::za4_gpu::tr70(fenv));
      const auto _interp178 = ntIm(DiFfRG::za4_gpu::tr71(fenv));
      const auto _interp180 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp184 = ntIm(DiFfRG::za4_gpu::tr72(fenv));
      const auto _interp185 = ntIm(DiFfRG::za4_gpu::tr59(fenv));
      const auto _interp187 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp191 = ntIm(DiFfRG::za4_gpu::tr60(fenv));
      const auto _interp192 = ntIm(DiFfRG::za4_gpu::tr57(fenv));
      const auto _interp193 = ntIm(DiFfRG::za4_gpu::tr58(fenv));
      const auto _interp194 = ntIm(DiFfRG::za4_gpu::tr61(fenv));
      const auto _interp196 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp197 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp204 = ntIm(DiFfRG::za4_gpu::tr62(fenv));
      const auto _interp205 = ntIm(DiFfRG::za4_gpu::tr75(fenv));
      const auto _interp207 = ntIm(DiFfRG::za4_gpu::tr76(fenv));
      const auto _interp208 = ntIm(DiFfRG::za4_gpu::tr77(fenv));
      const auto _den32 = powr<-1>(powr<2>(_interp174) + powr<2>(_interp162 * _interp175 + _interp176 * sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp214 = ntIm(DiFfRG::za4_gpu::tr78(fenv));
      const auto _interp215 = ntIm(DiFfRG::za4_gpu::tr73(fenv));
      const auto _interp216 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp217 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp218 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp219 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp220 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp221 = ntIm(DiFfRG::za4_gpu::tr74(fenv));
      const auto _interp222 = ntIm(DiFfRG::za4_gpu::tr79(fenv));
      const auto _interp223 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp224 = ntIm(DiFfRG::za4_gpu::tr80(fenv));
      const auto _interp225 = ntIm(DiFfRG::za4_gpu::tr63(fenv));
      const auto _interp226 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp227 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp228 = Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp229 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp230 = Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp231 = ntIm(DiFfRG::za4_gpu::tr64(fenv));
      const auto _interp232 = ntIm(DiFfRG::za4_gpu::tr65(fenv));
      const auto _interp233 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp234 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp235 = Mq(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp236 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp237 = Zq(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp238 = ntIm(DiFfRG::za4_gpu::tr66(fenv));
      const auto _interp239 = ntIm(DiFfRG::za4_gpu::tr67(fenv));
      const auto _interp240 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp241 = ntIm(DiFfRG::za4_gpu::tr68(fenv));
      const auto _den34 = powr<-1>(powr<2>(_interp218) + powr<2>(_interp162 * _interp219 + _interp220 * sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den35 = powr<-1>(powr<2>(_interp228) + powr<2>(_interp162 * _interp229 + _interp230 * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den36 = powr<-1>(powr<2>(_interp235) + powr<2>(_interp162 * _interp236 + _interp237 * sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      // clang-format off
      using _T = decltype(_den1 + _den10 + _den12 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den2 + _den20 + _den21 + _den22 + _den23 + _den24 + _den25 + _den26 + _den27 + _den28 + _den29 + _den3 + _den30 + _den31 + _den32 + _den33 + _den34 + _den35 + _den36 + _den37 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp122 + _interp123 + _interp124 + _interp125 + _interp126 + _interp127 + _interp128 + _interp129 + _interp13 + _interp130 + _interp131 + _interp132 + _interp133 + _interp134 + _interp135 + _interp136 + _interp137 + _interp138 + _interp139 + _interp14 + _interp140 + _interp141 + _interp142 + _interp143 + _interp144 + _interp145 + _interp146 + _interp147 + _interp148 + _interp149 + _interp15 + _interp150 + _interp151 + _interp152 + _interp153 + _interp154 + _interp155 + _interp156 + _interp157 + _interp158 + _interp159 + _interp16 + _interp160 + _interp161 + _interp162 + _interp163 + _interp164 + _interp165 + _interp166 + _interp167 + _interp168 + _interp169 + _interp17 + _interp170 + _interp171 + _interp172 + _interp173 + _interp174 + _interp175 + _interp176 + _interp177 + _interp178 + _interp179 + _interp18 + _interp180 + _interp181 + _interp182 + _interp183 + _interp184 + _interp185 + _interp186 + _interp187 + _interp188 + _interp189 + _interp19 + _interp190 + _interp191 + _interp192 + _interp193 + _interp194 + _interp195 + _interp196 + _interp197 + _interp198 + _interp199 + _interp2 + _interp20 + _interp200 + _interp201 + _interp202 + _interp203 + _interp204 + _interp205 + _interp206 + _interp207 + _interp208 + _interp209 + _interp21 + _interp210 + _interp211 + _interp212 + _interp213 + _interp214 + _interp215 + _interp216 + _interp217 + _interp218 + _interp219 + _interp22 + _interp220 + _interp221 + _interp222 + _interp223 + _interp224 + _interp225 + _interp226 + _interp227 + _interp228 + _interp229 + _interp23 + _interp230 + _interp231 + _interp232 + _interp233 + _interp234 + _interp235 + _interp236 + _interp237 + _interp238 + _interp239 + _interp24 + _interp240 + _interp241 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + cosl1p2 + k + l1 + p);
      // clang-format on
      _T _acc{};
      { // subkernel 1
        const auto _den11 = -powr<-1>(_interp100 * _interp41 + _interp124 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        // clang-format off
        _acc += fma(-0.0004251700680272108, _den10 * _den13 * _den21 * _den3 * _interp109 * _interp110 * _interp114 * _interp115 * _interp118 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(0.0004251700680272108, _den10 * _den13 * _den21 * _den3 * _interp109 * _interp110 * _interp114 * _interp115 * _interp119 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(-0.0004251700680272108, _den10 * _den13 * _den24 * _den3 * _interp109 * _interp114 * _interp143 * _interp144 * _interp145 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(0.0004251700680272108, _den10 * _den13 * _den24 * _den3 * _interp109 * _interp114 * _interp144 * _interp145 * _interp147 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4), fma(0.0004251700680272108, _den11 * _den22 * _den3 * _den9 * _interp120 * _interp121 * _interp122 * _interp123 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(-0.0004251700680272108, _den11 * _den22 * _den3 * _den9 * _interp121 * _interp122 * _interp123 * _interp126 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(-0.0004251700680272108, _den11 * _den25 * _den3 * _den9 * _interp121 * _interp148 * _interp149 * _interp150 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(0.0004251700680272108, _den11 * _den25 * _den3 * _den9 * _interp121 * _interp149 * _interp150 * _interp152 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(-0.0004251700680272108, _den13 * _den25 * _den3 * _den9 * _interp114 * _interp149 * _interp153 * _interp154 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(0.0004251700680272108, _den13 * _den25 * _den3 * _den9 * _interp114 * _interp149 * _interp154 * _interp155 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96, fma(0.0004251700680272108, _den12 * _den22 * _den26 * _den3 * _interp128 * _interp130 * _interp131 * _interp132 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(-0.0004251700680272108, _den12 * _den22 * _den26 * _den3 * _interp128 * _interp131 * _interp132 * _interp134 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(0.0004251700680272108, _den10 * _den23 * _den26 * _den3 * _interp109 * _interp135 * _interp136 * _interp137 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(-0.0004251700680272108, _den10 * _den23 * _den26 * _den3 * _interp109 * _interp136 * _interp137 * _interp139 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, 0.))))))))))))));
        // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(0.0004251700680272108, _den12 * _den23 * _den26 * _den3 * _interp131 * _interp136 * _interp140 * _interp141 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(-0.0004251700680272108, _den12 * _den23 * _den26 * _den3 * _interp131 * _interp136 * _interp141 * _interp142 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97, fma(0.0004251700680272108, _den22 * _den26 * _den3 * _den9 * _interp122 * _interp127 * _interp128 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp97, fma(-0.0004251700680272108, _den22 * _den26 * _den3 * _den9 * _interp122 * _interp128 * _interp129 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp97, fma(-0.0004251700680272108, _den13 * _den21 * _den3 * _den9 * _interp113 * _interp114 * _interp115 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp98, fma(0.0004251700680272108, _den13 * _den21 * _den3 * _den9 * _interp114 * _interp115 * _interp117 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp98, fma(0.0004251700680272108, _den10 * _den21 * _den26 * _den3 * _interp108 * _interp109 * _interp110 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97 * _interp99, fma(-0.0004251700680272108, _den10 * _den21 * _den26 * _den3 * _interp109 * _interp110 * _interp112 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp97 * _interp99, fma(-0.0004251700680272108, _den21 * _den26 * _den3 * _den9 * _interp107 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp96 * _interp97 * _interp98 * _interp99, fma(0.0004251700680272108, _den21 * _den26 * _den3 * _den9 * (_interp100 * _interp2 + (_interp101 + 50. * (-_interp100 + _interp102)) * _interp4) * _interp95 * _interp96 * _interp97 * _interp98 * _interp99, fma(0.0004251700680272108, _den17 * _den2 * _den20 * _den7 * _interp1 * _interp14 * _interp15 * _interp16 * _interp17 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den20 * _den7 * _interp14 * _interp15 * _interp18 * _interp21 * _interp22 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den20 * _den4 * _interp15 * _interp22 * _interp23 * _interp26 * _interp27 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den20 * _den5 * _interp15 * _interp16 * _interp28 * _interp31 * _interp32 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.0004251700680272108, _den15 * _den2 * _den20 * _den5 * _interp15 * _interp31 * _interp33 * _interp36 * _interp37 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den20 * _den4 * _interp15 * _interp26 * _interp37 * _interp38 * _interp39 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den4 * _den6 * _interp26 * _interp27 * _interp40 * _interp43 * _interp44 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den5 * _den8 * _interp31 * _interp36 * _interp45 * _interp48 * _interp49 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den4 * _den8 * _interp26 * _interp39 * _interp48 * _interp49 * _interp50 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * _den6 * _interp26 * _interp43 * _interp51 * _interp54 * _interp55 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den5 * _den8 * _interp31 * _interp48 * _interp56 * _interp59 * _interp60 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * _den8 * _interp26 * _interp48 * _interp54 * _interp61 * _interp62 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den5 * _interp31 * _interp36 * _interp63 * _interp64 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den4 * _interp26 * _interp39 * _interp64 * _interp65 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den20 * _interp15 * _interp37 * _interp64 * _interp66 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den8 * _interp48 * _interp49 * _interp64 * _interp67 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * powr<2>(_interp64) * _interp68 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den7 * _interp14 * _interp21 * _interp69 * _interp70 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(0.0004251700680272108, _den16 * _den2 * _den20 * _interp15 * _interp22 * _interp70 * _interp71 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den4 * _interp26 * _interp27 * _interp70 * _interp72 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * powr<2>(_interp70) * _interp73 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den20 * _interp15 * _interp16 * _interp74 * _interp75 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den7 * _interp14 * _interp17 * _interp75 * _interp76 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den5 * _interp31 * _interp32 * _interp75 * _interp77 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * powr<2>(_interp75) * _interp78 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den5 * _interp31 * _interp59 * _interp79 * _interp80 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * _interp26 * _interp54 * _interp81 * _interp82 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den6 * _interp26 * _interp43 * _interp83 * _interp84 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den20 * _den7 * _interp14 * _interp15 * _interp85 * _interp86 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den20 * _den5 * _interp15 * _interp31 * _interp87 * _interp88 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den20 * _den4 * _interp15 * _interp26 * _interp89 * _interp90 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den8 * _interp26 * _interp48 * _interp91 * _interp92 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 5
      // clang-format off
      _acc += fma(0.0004251700680272108, _den2 * _den5 * _den8 * _interp31 * _interp48 * _interp93 * _interp94 * (_interp2 * _interp3 + _interp4 * (_interp5 + 50. * _den1 * (-_interp3 + _interp6) * powr<6>(k))), fma(0.000850340136054422, _den14 * _den27 * _den32 * _den37 * _interp156 * _interp157 * _interp158 * _interp159 * _interp160 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den27 * _den32 * _den37 * _interp157 * _interp158 * _interp159 * _interp160 * _interp177 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den28 * _den32 * _den37 * _interp158 * _interp160 * _interp178 * _interp179 * _interp180 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den28 * _den32 * _den37 * _interp158 * _interp160 * _interp179 * _interp180 * _interp184 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den27 * _den31 * _den32 * _interp157 * _interp159 * _interp185 * _interp186 * _interp187 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den27 * _den31 * _den32 * _interp157 * _interp159 * _interp186 * _interp187 * _interp191 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den28 * _den31 * _den32 * _interp179 * _interp180 * _interp186 * _interp187 * _interp192 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den28 * _den31 * _den32 * _interp179 * _interp180 * _interp186 * _interp187 * _interp193 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den27 * _den29 * _den33 * _interp157 * _interp194 * _interp195 * _interp196 * _interp197 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den27 * _den29 * _den33 * _interp157 * _interp195 * _interp196 * _interp197 * _interp204 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den27 * _den33 * _den37 * _interp157 * _interp158 * _interp196 * _interp205 * _interp206 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den27 * _den33 * _den37 * _interp157 * _interp158 * _interp196 * _interp206 * _interp207 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den30 * _den33 * _den37 * _interp158 * _interp206 * _interp208 * _interp209 * _interp210 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 6
      // clang-format off
      _acc += fma(-0.000850340136054422, _den14 * _den30 * _den33 * _den37 * _interp158 * _interp206 * _interp209 * _interp210 * _interp214 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den28 * _den34 * _den37 * _interp158 * _interp179 * _interp215 * _interp216 * _interp217 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den28 * _den34 * _den37 * _interp158 * _interp179 * _interp216 * _interp217 * _interp221 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den30 * _den34 * _den37 * _interp158 * _interp209 * _interp216 * _interp222 * _interp223 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den30 * _den34 * _den37 * _interp158 * _interp209 * _interp216 * _interp223 * _interp224 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den28 * _den31 * _den35 * _interp179 * _interp186 * _interp225 * _interp226 * _interp227 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den28 * _den31 * _den35 * _interp179 * _interp186 * _interp226 * _interp227 * _interp231 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den27 * _den29 * _den36 * _interp157 * _interp195 * _interp232 * _interp233 * _interp234 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den27 * _den29 * _den36 * _interp157 * _interp195 * _interp233 * _interp234 * _interp238 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(-0.000850340136054422, _den14 * _den27 * _den31 * _den36 * _interp157 * _interp186 * _interp233 * _interp239 * _interp240 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), fma(0.000850340136054422, _den14 * _den27 * _den31 * _den36 * _interp157 * _interp186 * _interp233 * _interp240 * _interp241 * (-_interp161 * _interp162 * powr<-1>(l1) - _interp163 * (_interp164 + 50. * (-_interp162 + _interp165)) * powr<-1>(l1)), 0.)))))))))));
      // clang-format on
      }
      return _acc;
    }
    #else                              // the imaginary part survives: genuinely complex
    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      double fenv[(DiFfRG::za4_gpu::nenv) > 0 ? (DiFfRG::za4_gpu::nenv) : 1];
      const double dr_0 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_1 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_2 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_3 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_4 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_5 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_6 = Mq(l1);
      const double dr_7 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
      const double dr_8 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_9 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const double dr_10 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_11 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const double dr_12 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_13 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_14 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_15 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_16 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_17 = -sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_18 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_19 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_20 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_21 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p - 2. * cosl1p3 * l1 * p + 1.333333333333333 * powr<2>(p)));
      const double dr_22 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const double dr_23 = -sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))) * RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)) * Zq(k) - Zq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      DiFfRG::za4_gpu::fill(fenv, l1, cos1, cos2, phi, p, dr_0, dr_1, dr_2, dr_3, dr_4, dr_5, dr_6, dr_7, dr_8, dr_9, dr_10, dr_11, dr_12, dr_13, dr_14, dr_15, dr_16, dr_17, dr_18, dr_19, dr_20, dr_21, dr_22, dr_23);
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
      const auto _interp52 = ZA4(0.7071067811865475 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
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
      const auto _interp99 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp100 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp101 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1))))));
      const auto _interp102 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-1., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp103 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _interp104 = Zq(k);
      const auto _interp105 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp106 = dtZq(k);
      const auto _interp107 = Zq(1.02 * k);
      const auto _interp108 = Mq(l1);
      const auto _interp109 = Zq(l1);
      const auto _interp110 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp111 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp112 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp113 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp114 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)))));
      const auto _interp115 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))))));
      const auto _interp116 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp117 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp118 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p2 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp124 = ZAqbq1(0.816496580927726 * sqrt(fma(cosl1p3, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp125 = ZAqbq1(0.816496580927726 * sqrt(fma(-1., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(cosl1p3, l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
      const auto _interp126 = Mq(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp127 = RF(powr<2>(k), fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp128 = Zq(sqrt(fma(2., cosl1p3 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp135 = Mq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp136 = RF(powr<2>(k), fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1)))));
      const auto _interp137 = Zq(sqrt(fma(-2., cosl1p1 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.333333333333333, powr<2>(p), powr<2>(l1))))));
      const auto _interp138 = ZAqbq1(0.816496580927726 * sqrt(fma(-2., cosl1p1 * l1 * p, fma(-1., cosl1p2 * l1 * p, fma(-2., cosl1p3 * l1 * p, fma(1.666666666666667, powr<2>(p), powr<2>(l1)))))));
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
      const auto _den14 = powr<-2>(fma(powr<-2>(l1) * powr<2>(_interp104 * _interp105 + _interp109 * l1), powr<2>(l1), powr<2>(_interp108)));
      const auto _den16 = powr<-1>(fma(_interp17, _interp2, fma(_interp18, powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den17 = powr<-1>(fma(_interp11, _interp2, fma(_interp12, powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den18 = powr<-1>(fma(_interp2, _interp46, fma(_interp47, powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den19 = powr<-1>(fma(_interp2, _interp42, fma(_interp43, powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), 0.)));
      const auto _den20 = powr<-1>(fma(_interp2, _interp9, fma(_interp10, powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), 0.)));
      const auto _den26 = powr<-1>(fma(-1., _interp66 * _interp9, fma(-1., _interp71 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)), 0.)));
      const auto _den27 = powr<-1>(fma(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p), powr<2>(_interp112 + _interp104 * _interp111 * sqrt(powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))), powr<2>(_interp110)));
      const auto _den31 = powr<-1>(fma(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p), powr<2>(_interp128 + _interp104 * _interp127 * sqrt(powr<-1>(powr<2>(l1) + 2. * cosl1p3 * l1 * p + powr<2>(p)))), powr<2>(_interp126)));
      const auto _den32 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp118 + _interp104 * _interp117 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp116)));
      const auto _den33 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p), powr<2>(_interp137 + _interp104 * _interp136 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)))), powr<2>(_interp135)));
      const auto _den37 = powr<-1>(fma(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p), powr<2>(_interp115 + _interp104 * _interp114 * sqrt(powr<-1>(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)))), powr<2>(_interp113)));
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
      const auto _den21 = -3. * powr<-1>(3. * _interp29 * _interp66 + _interp72 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den22 = -3. * powr<-1>(3. * _interp17 * _interp66 + _interp83 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _interp88 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp89 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp90 = Zc(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp91 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp92 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp93 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp94 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp95 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp96 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp97 = Zc(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp98 = ZAcbc(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _den23 = -3. * powr<-1>(3. * _interp11 * _interp66 + _interp90 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den24 = -3. * powr<-1>(3. * _interp46 * _interp66 + _interp94 * (3. * powr<2>(l1) + 6. * (cosl1p1 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den25 = -3. * powr<-1>(3. * _interp42 * _interp66 + _interp97 * (3. * powr<2>(l1) + 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _interp30 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp31 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp32 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp33 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp37 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp41 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2 - cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp45 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp49 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp50 = ZA3(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp51 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _den15 = 3. * powr<-1>(3. * _interp2 * _interp29 + _interp30 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _interp53 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp54 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p1 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp55 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) + 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp56 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p2) * l1 * p + 2. * powr<2>(p)));
      const auto _interp57 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + 2. * cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp58 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp59 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (2. * cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp60 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p1 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp61 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) + (-cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp119 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp120 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp121 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp122 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp123 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp139 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp140 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp141 = Mq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp142 = RF(powr<2>(k), powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p));
      const auto _interp143 = Zq(sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p)));
      const auto _interp144 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 + 2. * (cosl1p2 + cosl1p3)) * l1 * p + 10. * powr<2>(p)));
      const auto _interp145 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp146 = Mq(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp147 = RF(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp148 = Zq(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp149 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp150 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp151 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p1 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp152 = Mq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp153 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp154 = Zq(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _den28 = powr<-1>(powr<2>(_interp121) + powr<2>(_interp104 * _interp122 + _interp123 * sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den30 = powr<-1>(powr<2>(_interp141) + powr<2>(_interp104 * _interp142 + _interp143 * sqrt(powr<2>(l1) - 2. * cosl1p3 * l1 * p + powr<2>(p))));
      const auto _den34 = powr<-1>(powr<2>(_interp146) + powr<2>(_interp104 * _interp147 + _interp148 * sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _den35 = powr<-1>(powr<2>(_interp152) + powr<2>(_interp104 * _interp153 + _interp154 * sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      const auto _interp129 = ZAqbq1(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp130 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp131 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (cosl1p1 - cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp132 = Mq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp133 = RF(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp134 = Zq(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp155 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (-cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp156 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (2. * cosl1p2 + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp157 = Mq(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp158 = RF(powr<2>(k), powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp159 = Zq(sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp160 = ZAqbq1(0.3333333333333333 * sqrt(6. * powr<2>(l1) + 6. * (cosl1p2 + 2. * cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _den29 = powr<-1>(powr<2>(_interp132) + powr<2>(_interp104 * _interp133 + _interp134 * sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p))));
      const auto _den36 = powr<-1>(powr<2>(_interp157) + powr<2>(_interp104 * _interp158 + _interp159 * sqrt(powr<2>(l1) + 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p))));
      // clang-format off
      using _T = decltype(complex<double>(0.,0.000850340136054422) + _den1 + _den10 + _den11 + _den12 + _den13 + _den14 + _den15 + _den16 + _den17 + _den18 + _den19 + _den2 + _den20 + _den21 + _den22 + _den23 + _den24 + _den25 + _den26 + _den27 + _den28 + _den29 + _den3 + _den30 + _den31 + _den32 + _den33 + _den34 + _den35 + _den36 + _den37 + _den4 + _den5 + _den6 + _den7 + _den8 + _den9 + DiFfRG::za4_gpu::tr0(fenv) + DiFfRG::za4_gpu::tr10(fenv) + DiFfRG::za4_gpu::tr11(fenv) + DiFfRG::za4_gpu::tr12(fenv) + DiFfRG::za4_gpu::tr13(fenv) + DiFfRG::za4_gpu::tr14(fenv) + DiFfRG::za4_gpu::tr15(fenv) + DiFfRG::za4_gpu::tr16(fenv) + DiFfRG::za4_gpu::tr17(fenv) + DiFfRG::za4_gpu::tr18(fenv) + DiFfRG::za4_gpu::tr19(fenv) + DiFfRG::za4_gpu::tr1(fenv) + DiFfRG::za4_gpu::tr20(fenv) + DiFfRG::za4_gpu::tr21(fenv) + DiFfRG::za4_gpu::tr22(fenv) + DiFfRG::za4_gpu::tr23(fenv) + DiFfRG::za4_gpu::tr24(fenv) + DiFfRG::za4_gpu::tr25(fenv) + DiFfRG::za4_gpu::tr26(fenv) + DiFfRG::za4_gpu::tr27(fenv) + DiFfRG::za4_gpu::tr28(fenv) + DiFfRG::za4_gpu::tr29(fenv) + DiFfRG::za4_gpu::tr2(fenv) + DiFfRG::za4_gpu::tr30(fenv) + DiFfRG::za4_gpu::tr31(fenv) + DiFfRG::za4_gpu::tr32(fenv) + DiFfRG::za4_gpu::tr33(fenv) + DiFfRG::za4_gpu::tr34(fenv) + DiFfRG::za4_gpu::tr35(fenv) + DiFfRG::za4_gpu::tr36(fenv) + DiFfRG::za4_gpu::tr37(fenv) + DiFfRG::za4_gpu::tr38(fenv) + DiFfRG::za4_gpu::tr39(fenv) + DiFfRG::za4_gpu::tr3(fenv) + DiFfRG::za4_gpu::tr40(fenv) + DiFfRG::za4_gpu::tr41(fenv) + DiFfRG::za4_gpu::tr42(fenv) + DiFfRG::za4_gpu::tr43(fenv) + DiFfRG::za4_gpu::tr44(fenv) + DiFfRG::za4_gpu::tr45(fenv) + DiFfRG::za4_gpu::tr46(fenv) + DiFfRG::za4_gpu::tr47(fenv) + DiFfRG::za4_gpu::tr48(fenv) + DiFfRG::za4_gpu::tr49(fenv) + DiFfRG::za4_gpu::tr4(fenv) + DiFfRG::za4_gpu::tr50(fenv) + DiFfRG::za4_gpu::tr51(fenv) + DiFfRG::za4_gpu::tr52(fenv) + DiFfRG::za4_gpu::tr53(fenv) + DiFfRG::za4_gpu::tr54(fenv) + DiFfRG::za4_gpu::tr55(fenv) + DiFfRG::za4_gpu::tr56(fenv) + DiFfRG::za4_gpu::tr57(fenv) + DiFfRG::za4_gpu::tr58(fenv) + DiFfRG::za4_gpu::tr59(fenv) + DiFfRG::za4_gpu::tr5(fenv) + DiFfRG::za4_gpu::tr60(fenv) + DiFfRG::za4_gpu::tr61(fenv) + DiFfRG::za4_gpu::tr62(fenv) + DiFfRG::za4_gpu::tr63(fenv) + DiFfRG::za4_gpu::tr64(fenv) + DiFfRG::za4_gpu::tr65(fenv) + DiFfRG::za4_gpu::tr66(fenv) + DiFfRG::za4_gpu::tr67(fenv) + DiFfRG::za4_gpu::tr68(fenv) + DiFfRG::za4_gpu::tr69(fenv) + DiFfRG::za4_gpu::tr6(fenv) + DiFfRG::za4_gpu::tr70(fenv) + DiFfRG::za4_gpu::tr71(fenv) + DiFfRG::za4_gpu::tr72(fenv) + DiFfRG::za4_gpu::tr73(fenv) + DiFfRG::za4_gpu::tr74(fenv) + DiFfRG::za4_gpu::tr75(fenv) + DiFfRG::za4_gpu::tr76(fenv) + DiFfRG::za4_gpu::tr77(fenv) + DiFfRG::za4_gpu::tr78(fenv) + DiFfRG::za4_gpu::tr79(fenv) + DiFfRG::za4_gpu::tr7(fenv) + DiFfRG::za4_gpu::tr80(fenv) + DiFfRG::za4_gpu::tr8(fenv) + DiFfRG::za4_gpu::tr9(fenv) + _interp1 + _interp10 + _interp100 + _interp101 + _interp102 + _interp103 + _interp104 + _interp105 + _interp106 + _interp107 + _interp108 + _interp109 + _interp11 + _interp110 + _interp111 + _interp112 + _interp113 + _interp114 + _interp115 + _interp116 + _interp117 + _interp118 + _interp119 + _interp12 + _interp120 + _interp121 + _interp122 + _interp123 + _interp124 + _interp125 + _interp126 + _interp127 + _interp128 + _interp129 + _interp13 + _interp130 + _interp131 + _interp132 + _interp133 + _interp134 + _interp135 + _interp136 + _interp137 + _interp138 + _interp139 + _interp14 + _interp140 + _interp141 + _interp142 + _interp143 + _interp144 + _interp145 + _interp146 + _interp147 + _interp148 + _interp149 + _interp15 + _interp150 + _interp151 + _interp152 + _interp153 + _interp154 + _interp155 + _interp156 + _interp157 + _interp158 + _interp159 + _interp16 + _interp160 + _interp17 + _interp18 + _interp19 + _interp2 + _interp20 + _interp21 + _interp22 + _interp23 + _interp24 + _interp25 + _interp26 + _interp27 + _interp28 + _interp29 + _interp3 + _interp30 + _interp31 + _interp32 + _interp33 + _interp34 + _interp35 + _interp36 + _interp37 + _interp38 + _interp39 + _interp4 + _interp40 + _interp41 + _interp42 + _interp43 + _interp44 + _interp45 + _interp46 + _interp47 + _interp48 + _interp49 + _interp5 + _interp50 + _interp51 + _interp52 + _interp53 + _interp54 + _interp55 + _interp56 + _interp57 + _interp58 + _interp59 + _interp6 + _interp60 + _interp61 + _interp62 + _interp63 + _interp64 + _interp65 + _interp66 + _interp67 + _interp68 + _interp69 + _interp7 + _interp70 + _interp71 + _interp72 + _interp73 + _interp74 + _interp75 + _interp76 + _interp77 + _interp78 + _interp79 + _interp8 + _interp80 + _interp81 + _interp82 + _interp83 + _interp84 + _interp85 + _interp86 + _interp87 + _interp88 + _interp89 + _interp9 + _interp90 + _interp91 + _interp92 + _interp93 + _interp94 + _interp95 + _interp96 + _interp97 + _interp98 + _interp99 + k + l1);
      // clang-format on
      _T _acc{};
      { // subkernel 1
      // clang-format off
      _acc += fma(0.0004251700680272108, _den21 * _den26 * _den3 * _den9 * DiFfRG::za4_gpu::tr45(fenv) * _interp62 * _interp63 * _interp64 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))), fma(-0.0004251700680272108, _den21 * _den26 * _den3 * _den9 * DiFfRG::za4_gpu::tr46(fenv) * _interp62 * _interp63 * _interp64 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))), fma(0.0004251700680272108, _den10 * _den21 * _den26 * _den3 * DiFfRG::za4_gpu::tr47(fenv) * _interp63 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74, fma(-0.0004251700680272108, _den10 * _den21 * _den26 * _den3 * DiFfRG::za4_gpu::tr48(fenv) * _interp63 * _interp65 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74, fma(-0.0004251700680272108, _den13 * _den21 * _den3 * _den9 * DiFfRG::za4_gpu::tr33(fenv) * _interp62 * _interp64 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp77, fma(0.0004251700680272108, _den13 * _den21 * _den3 * _den9 * DiFfRG::za4_gpu::tr34(fenv) * _interp62 * _interp64 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp77, fma(-0.0004251700680272108, _den10 * _den13 * _den21 * _den3 * DiFfRG::za4_gpu::tr35(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74 * _interp76 * _interp77, fma(0.0004251700680272108, _den10 * _den13 * _den21 * _den3 * DiFfRG::za4_gpu::tr36(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp74 * _interp76 * _interp77, fma(0.0004251700680272108, _den11 * _den22 * _den3 * _den9 * DiFfRG::za4_gpu::tr37(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp80 * _interp81, fma(-0.0004251700680272108, _den11 * _den22 * _den3 * _den9 * DiFfRG::za4_gpu::tr38(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp80 * _interp81, fma(0.0004251700680272108, _den22 * _den26 * _den3 * _den9 * DiFfRG::za4_gpu::tr49(fenv) * _interp62 * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp80 * _interp84, fma(-0.0004251700680272108, _den22 * _den26 * _den3 * _den9 * DiFfRG::za4_gpu::tr50(fenv) * _interp62 * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp80 * _interp84, fma(0.0004251700680272108, _den12 * _den22 * _den26 * _den3 * DiFfRG::za4_gpu::tr51(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp84 * _interp85 * _interp86, fma(-0.0004251700680272108, _den12 * _den22 * _den26 * _den3 * DiFfRG::za4_gpu::tr52(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp84 * _interp85 * _interp86, 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 2
      // clang-format off
      _acc += fma(0.0004251700680272108, _den10 * _den23 * _den26 * _den3 * DiFfRG::za4_gpu::tr53(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp88 * _interp89, fma(-0.0004251700680272108, _den10 * _den23 * _den26 * _den3 * DiFfRG::za4_gpu::tr54(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp88 * _interp89, fma(0.0004251700680272108, _den12 * _den23 * _den26 * _den3 * DiFfRG::za4_gpu::tr55(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp85 * _interp88 * _interp91, fma(-0.0004251700680272108, _den12 * _den23 * _den26 * _den3 * DiFfRG::za4_gpu::tr56(fenv) * _interp63 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp85 * _interp88 * _interp91, fma(-0.0004251700680272108, _den10 * _den13 * _den24 * _den3 * DiFfRG::za4_gpu::tr39(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp76 * _interp92 * _interp93, fma(0.0004251700680272108, _den10 * _den13 * _den24 * _den3 * DiFfRG::za4_gpu::tr40(fenv) * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp73 * _interp76 * _interp92 * _interp93, fma(-0.0004251700680272108, _den11 * _den25 * _den3 * _den9 * DiFfRG::za4_gpu::tr41(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp95 * _interp96, fma(0.0004251700680272108, _den11 * _den25 * _den3 * _den9 * DiFfRG::za4_gpu::tr42(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp79 * _interp95 * _interp96, fma(-0.0004251700680272108, _den13 * _den25 * _den3 * _den9 * DiFfRG::za4_gpu::tr43(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp95 * _interp98, fma(0.0004251700680272108, _den13 * _den25 * _den3 * _den9 * DiFfRG::za4_gpu::tr44(fenv) * _interp62 * (_interp1 * _interp66 + _interp3 * (_interp67 + 50. * (-_interp66 + _interp68))) * _interp76 * _interp95 * _interp98, fma(0.0004251700680272108, _den17 * _den2 * _den20 * _den7 * DiFfRG::za4_gpu::tr10(fenv) * _interp13 * _interp14 * _interp15 * _interp16 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den20 * _den7 * DiFfRG::za4_gpu::tr8(fenv) * _interp13 * _interp14 * _interp19 * _interp20 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den20 * _den4 * DiFfRG::za4_gpu::tr9(fenv) * _interp14 * _interp20 * _interp23 * _interp24 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den20 * _den5 * DiFfRG::za4_gpu::tr11(fenv) * _interp14 * _interp15 * _interp27 * _interp28 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 3
      // clang-format off
      _acc += fma(0.0004251700680272108, _den15 * _den2 * _den20 * _den5 * DiFfRG::za4_gpu::tr6(fenv) * _interp14 * _interp27 * _interp31 * _interp32 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den20 * _den4 * DiFfRG::za4_gpu::tr7(fenv) * _interp14 * _interp23 * _interp32 * _interp33 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den4 * _den6 * DiFfRG::za4_gpu::tr2(fenv) * _interp23 * _interp24 * _interp36 * _interp37 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den5 * _den8 * DiFfRG::za4_gpu::tr0(fenv) * _interp27 * _interp31 * _interp40 * _interp41 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den4 * _den8 * DiFfRG::za4_gpu::tr1(fenv) * _interp23 * _interp33 * _interp40 * _interp41 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * _den6 * DiFfRG::za4_gpu::tr4(fenv) * _interp23 * _interp36 * _interp44 * _interp45 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den5 * _den8 * DiFfRG::za4_gpu::tr3(fenv) * _interp27 * _interp40 * _interp48 * _interp49 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * _den8 * DiFfRG::za4_gpu::tr5(fenv) * _interp23 * _interp40 * _interp44 * _interp50 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den5 * DiFfRG::za4_gpu::tr24(fenv) * _interp27 * _interp31 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den20 * DiFfRG::za4_gpu::tr13(fenv) * _interp14 * _interp32 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den4 * DiFfRG::za4_gpu::tr25(fenv) * _interp23 * _interp33 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * _den8 * DiFfRG::za4_gpu::tr12(fenv) * _interp40 * _interp41 * _interp51 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den15 * _den2 * DiFfRG::za4_gpu::tr26(fenv) * powr<2>(_interp51) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den7 * DiFfRG::za4_gpu::tr27(fenv) * _interp13 * _interp19 * _interp52 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 4
      // clang-format off
      _acc += fma(0.0004251700680272108, _den16 * _den2 * _den20 * DiFfRG::za4_gpu::tr14(fenv) * _interp14 * _interp20 * _interp52 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * _den4 * DiFfRG::za4_gpu::tr28(fenv) * _interp23 * _interp24 * _interp52 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den16 * _den2 * DiFfRG::za4_gpu::tr29(fenv) * powr<2>(_interp52) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den20 * DiFfRG::za4_gpu::tr16(fenv) * _interp14 * _interp15 * _interp53 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den7 * DiFfRG::za4_gpu::tr30(fenv) * _interp13 * _interp16 * _interp53 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * _den5 * DiFfRG::za4_gpu::tr31(fenv) * _interp27 * _interp28 * _interp53 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den17 * _den2 * DiFfRG::za4_gpu::tr32(fenv) * powr<2>(_interp53) * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den18 * _den2 * _den5 * DiFfRG::za4_gpu::tr15(fenv) * _interp27 * _interp48 * _interp54 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den19 * _den2 * _den4 * DiFfRG::za4_gpu::tr17(fenv) * _interp23 * _interp44 * _interp55 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den6 * DiFfRG::za4_gpu::tr18(fenv) * _interp23 * _interp36 * _interp56 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den20 * _den7 * DiFfRG::za4_gpu::tr19(fenv) * _interp13 * _interp14 * _interp57 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den20 * _den5 * DiFfRG::za4_gpu::tr20(fenv) * _interp14 * _interp27 * _interp58 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den20 * _den4 * DiFfRG::za4_gpu::tr21(fenv) * _interp14 * _interp23 * _interp59 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(0.0004251700680272108, _den2 * _den4 * _den8 * DiFfRG::za4_gpu::tr22(fenv) * _interp23 * _interp40 * _interp60 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 5
      // clang-format off
      _acc += fma(0.0004251700680272108, _den2 * _den5 * _den8 * DiFfRG::za4_gpu::tr23(fenv) * _interp27 * _interp40 * _interp61 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den28 * _den32 * _den37 * DiFfRG::za4_gpu::tr71(fenv) * _interp100 * _interp102 * _interp119 * _interp120 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den28 * _den32 * _den37 * DiFfRG::za4_gpu::tr72(fenv) * _interp100 * _interp102 * _interp119 * _interp120 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den28 * _den31 * _den32 * DiFfRG::za4_gpu::tr57(fenv) * _interp119 * _interp120 * _interp124 * _interp125 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den28 * _den31 * _den32 * DiFfRG::za4_gpu::tr58(fenv) * _interp119 * _interp120 * _interp124 * _interp125 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den30 * _den33 * _den37 * DiFfRG::za4_gpu::tr77(fenv) * _interp100 * _interp138 * _interp139 * _interp140 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den30 * _den33 * _den37 * DiFfRG::za4_gpu::tr78(fenv) * _interp100 * _interp138 * _interp139 * _interp140 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den28 * _den34 * _den37 * DiFfRG::za4_gpu::tr73(fenv) * _interp100 * _interp119 * _interp144 * _interp145 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den28 * _den34 * _den37 * DiFfRG::za4_gpu::tr74(fenv) * _interp100 * _interp119 * _interp144 * _interp145 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den30 * _den34 * _den37 * DiFfRG::za4_gpu::tr79(fenv) * _interp100 * _interp139 * _interp144 * _interp149 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den30 * _den34 * _den37 * DiFfRG::za4_gpu::tr80(fenv) * _interp100 * _interp139 * _interp144 * _interp149 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den28 * _den31 * _den35 * DiFfRG::za4_gpu::tr63(fenv) * _interp119 * _interp124 * _interp150 * _interp151 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den28 * _den31 * _den35 * DiFfRG::za4_gpu::tr64(fenv) * _interp119 * _interp124 * _interp150 * _interp151 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den27 * _den32 * _den37 * DiFfRG::za4_gpu::tr69(fenv) * _interp100 * _interp101 * _interp102 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), 0.))))))))))))));
      // clang-format on
      }
      { // subkernel 6
      // clang-format off
      _acc += fma(complex<double>(0.,0.000850340136054422), _den14 * _den27 * _den32 * _den37 * DiFfRG::za4_gpu::tr70(fenv) * _interp100 * _interp101 * _interp102 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den27 * _den31 * _den32 * DiFfRG::za4_gpu::tr59(fenv) * _interp101 * _interp124 * _interp125 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den27 * _den31 * _den32 * DiFfRG::za4_gpu::tr60(fenv) * _interp101 * _interp124 * _interp125 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den27 * _den29 * _den33 * DiFfRG::za4_gpu::tr61(fenv) * _interp129 * _interp130 * _interp131 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den27 * _den29 * _den33 * DiFfRG::za4_gpu::tr62(fenv) * _interp129 * _interp130 * _interp131 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den27 * _den33 * _den37 * DiFfRG::za4_gpu::tr75(fenv) * _interp100 * _interp130 * _interp138 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den27 * _den33 * _den37 * DiFfRG::za4_gpu::tr76(fenv) * _interp100 * _interp130 * _interp138 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den27 * _den29 * _den36 * DiFfRG::za4_gpu::tr65(fenv) * _interp129 * _interp155 * _interp156 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den27 * _den29 * _den36 * DiFfRG::za4_gpu::tr66(fenv) * _interp129 * _interp155 * _interp156 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,0.000850340136054422), _den14 * _den27 * _den31 * _den36 * DiFfRG::za4_gpu::tr67(fenv) * _interp124 * _interp155 * _interp160 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), fma(complex<double>(0.,-0.000850340136054422), _den14 * _den27 * _den31 * _den36 * DiFfRG::za4_gpu::tr68(fenv) * _interp124 * _interp155 * _interp160 * _interp99 * (-_interp103 * _interp104 * powr<-1>(l1) - _interp105 * (_interp106 + 50. * (-_interp104 + _interp107)) * powr<-1>(l1)), 0.)))))))))));
      // clang-format on
      }
      return _acc;
    }
    #endif

    // clang-format off
    static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Mq)
    // clang-format on
    {
      return 0.;
    }
    private:
    static KOKKOS_INLINE_FUNCTION double ntRe(double x) { return x; }
    template <class T> static KOKKOS_INLINE_FUNCTION double ntRe(const T &z) { return z.real(); }
    static KOKKOS_INLINE_FUNCTION double ntIm(double) { return 0.0; }
    template <class T> static KOKKOS_INLINE_FUNCTION double ntIm(const T &z) { return z.imag(); }
  };
}
using DiFfRG::ZA4_gpu_kernel;