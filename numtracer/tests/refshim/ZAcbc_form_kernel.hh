#pragma once

#include "shim.hpp"

namespace DiFfRG {
  template<typename _Regulator>
  class ZAcbc_form_kernel
  {
    public:
    using Regulator = _Regulator;

    static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;
      const double cosl1p1 = cos1;
      const double cosl1p2 = 0.5 * (-cos1 + sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      const double cosl1p3 = 0.5 * (-cos1 - sqrt(3. - 3. * powr<2>(cos1)) * cos2);
      // clang-format off
using _T = decltype(-0.5 * p * powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) * (4. * cosl1p2 * powr<3>(l1) - 8. * powr<3>(cosl1p2) * powr<3>(l1) + 6. * powr<2>(l1) * p - 4. * powr<3>(cosl1p1) * cosl1p2 * powr<2>(l1) * p - 4. * powr<2>(cosl1p2) * powr<2>(l1) * p - 2. * cosl1p2 * l1 * powr<2>(p) + 3. * powr<3>(p) - 2. * powr<2>(cosl1p1) * l1 * (2. * cosl1p2 * powr<2>(l1) + l1 * p + 6. * powr<2>(cosl1p2) * l1 * p - 5. * cosl1p2 * powr<2>(p)) + cosl1p1 * ((2. - 12. * powr<2>(cosl1p2)) * powr<3>(l1) + (5. - 4. * powr<2>(cosl1p2)) * l1 * powr<2>(p) + 6. * cosl1p2 * powr<3>(p) + powr<2>(l1) * (6. * cosl1p2 * p - 8. * powr<3>(cosl1p2) * p))) * (RBdot(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + RB(powr<2>(k), powr<2>(l1)) * (dtZA(pow(1. + powr<6>(k),0.16666666666666666667)) + 50. * powr<6>(k) * powr<-1>(1. + powr<6>(k)) * (-ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667))))) * powr<-2>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<2>(l1) * ZA(l1)) * powr<-1>(RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + (powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) * ZA(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)))) * ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p1 * l1 * p + powr<2>(p))) * ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 - cosl1p2) * l1 * p + powr<2>(p))) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p))) * powr<-1>(RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))) + 0.5 * powr<-1>(1. + powr<6>(k)) * p * powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * (-4. * cosl1p2 * powr<3>(l1) + 8. * powr<3>(cosl1p2) * powr<3>(l1) - 6. * powr<2>(l1) * p + 4. * powr<2>(cosl1p2) * powr<2>(l1) * p + powr<3>(cosl1p1) * l1 * (-4. * cosl1p2 * l1 - 14. * p) * p + 2. * cosl1p2 * l1 * powr<2>(p) - 3. * powr<3>(p) + powr<2>(cosl1p1) * (4. * cosl1p2 * powr<3>(l1) + 12. * powr<2>(l1) * p - 12. * powr<2>(cosl1p2) * powr<2>(l1) * p - 18. * cosl1p2 * l1 * powr<2>(p) + 6. * powr<3>(p)) + cosl1p1 * ((-2. + 12. * powr<2>(cosl1p2)) * powr<3>(l1) + cosl1p2 * (14. - 8. * powr<2>(cosl1p2)) * powr<2>(l1) * p + (7. - 4. * powr<2>(cosl1p2)) * l1 * powr<2>(p) + 6. * cosl1p2 * powr<3>(p))) * ((1. + 1. * powr<6>(k)) * dtZA(pow(1. + powr<6>(k),0.16666666666666666667)) * RB(powr<2>(k), powr<2>(l1)) + (1. + powr<6>(k)) * RBdot(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<6>(k) * RB(powr<2>(k), powr<2>(l1)) * (-50. * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + 50. * ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667)))) * powr<-2>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<2>(l1) * ZA(l1)) * powr<-1>(RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)) * ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)))) * ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p))) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * powr<-1>(RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) - 0.25 * (-1. + 2. * cosl1p1 * cosl1p2 + 2. * powr<2>(cosl1p2)) * powr<-1>(1. + powr<6>(k)) * (2. * cosl1p1 * l1 + 4. * cosl1p2 * l1 - 3. * p) * p * ((1. + 1. * powr<6>(k)) * dtZA(pow(1. + powr<6>(k),0.16666666666666666667)) * RB(powr<2>(k), powr<2>(l1)) + (1. + powr<6>(k)) * RBdot(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<6>(k) * RB(powr<2>(k), powr<2>(l1)) * (-50. * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + 50. * ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667)))) * powr<-2>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<2>(l1) * ZA(l1)) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p))) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p))) * powr<-1>(RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)))) * powr<-1>(RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) - 2. * cosl1p1 * l1 * p - 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) + 3. * p * powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) * (0.6666666666666666 * cosl1p2 * powr<3>(l1) - 1.333333333333333 * powr<3>(cosl1p2) * powr<3>(l1) - powr<2>(l1) * p + 0.6666666666666666 * powr<2>(cosl1p2) * powr<2>(l1) * p - 0.3333333333333333 * cosl1p2 * l1 * powr<2>(p) - 0.5 * powr<3>(p) + powr<3>(cosl1p1) * l1 * p * (-0.6666666666666666 * cosl1p2 * l1 + 2.333333333333333 * p) + powr<2>(cosl1p1) * (-0.6666666666666666 * cosl1p2 * powr<3>(l1) + 2. * powr<2>(l1) * p - 2. * powr<2>(cosl1p2) * powr<2>(l1) * p + 3. * cosl1p2 * l1 * powr<2>(p) + 1. * powr<3>(p)) + cosl1p1 * ((0.3333333333333333 - 2. * powr<2>(cosl1p2)) * powr<3>(l1) + cosl1p2 * (2.333333333333333 - 1.333333333333333 * powr<2>(cosl1p2)) * powr<2>(l1) * p + (-1.166666666666667 + 0.6666666666666666 * powr<2>(cosl1p2)) * l1 * powr<2>(p) + 1. * cosl1p2 * powr<3>(p))) * powr<-1>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<2>(l1) * ZA(l1)) * powr<-1>(RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + (powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)) * ZA(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)))) * ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p1 * l1 * p + powr<2>(p))) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * (dtZc(k) * RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + RBdot(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * (-50. * Zc(k) + 50. * Zc(1.02 * k))) * powr<-2>(RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) + 1. * (-0.5 + 1. * cosl1p1 * cosl1p2 + 1. * powr<2>(cosl1p2)) * p * (1. * cosl1p1 * l1 + 2. * cosl1p2 * l1 + 1.5 * p) * powr<-1>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<2>(l1) * ZA(l1)) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p))) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p))) * (dtZc(k) * RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) + RBdot(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * (-50. * Zc(k) + 50. * Zc(1.02 * k))) * powr<-1>(RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) * powr<-2>(RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) + 2. * cosl1p1 * l1 * p + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))) + 0.25 * (-1. + 2. * cosl1p1 * cosl1p2 + 2. * powr<2>(cosl1p2)) * p * (2. * cosl1p1 * l1 + 4. * cosl1p2 * l1 + 3. * p) * powr<-1>(RB(powr<2>(k), powr<2>(l1)) * ZA(pow(1. + powr<6>(k),0.16666666666666666667)) + powr<2>(l1) * ZA(l1)) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p))) * ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p))) * ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p))) * (RBdot(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(k) + RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * (dtZc(k) - 50. * Zc(k) + 50. * Zc(1.02 * k))) * powr<-2>(RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)))) * powr<-1>(RB(powr<2>(k), powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(k) + (powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)) * Zc(sqrt(powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)))));
      // clang-format on
      const auto _interp1 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _interp2 = ZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp3 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp4 = dtZA(pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp5 = ZA(1.02 * pow(1. + powr<6>(k),0.16666666666666666667));
      const auto _interp6 = ZA(l1);
      const auto _interp11 = ZAcbc(0.816496580927726 * sqrt(fma(-1., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp12 = RB(powr<2>(k), fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp13 = Zc(k);
      const auto _interp14 = Zc(sqrt(fma(-2., cosl1p2 * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp18 = ZAcbc(0.816496580927726 * sqrt(fma(-1., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp20 = RB(powr<2>(k), fma(-2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp21 = Zc(sqrt(fma(-2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp23 = ZAcbc(0.816496580927726 * sqrt(fma(cosl1p1 + cosl1p2, l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _interp25 = dtZc(k);
      const auto _interp26 = RB(powr<2>(k), fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp27 = RBdot(powr<2>(k), fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p)));
      const auto _interp28 = Zc(1.02 * k);
      const auto _interp29 = Zc(sqrt(fma(2., (cosl1p1 + cosl1p2) * l1 * p, powr<2>(l1) + powr<2>(p))));
      const auto _den2 = powr<-1>(1. + powr<6>(k));
      const auto _den3 = powr<-2>(fma(_interp2, _interp3, fma(_interp6, powr<2>(l1), 0.)));
      const auto _den4 = powr<-1>(fma(_interp2, _interp3, fma(_interp6, powr<2>(l1), 0.)));
      const auto _den16 = powr<-1>(fma(_interp13, _interp20, fma(_interp21, powr<2>(l1) + (-2. * cosl1p1 - 2. * cosl1p2) * l1 * p + powr<2>(p), 0.)));
      const auto _den17 = powr<-2>(fma(_interp13, _interp26, fma(_interp29, powr<2>(l1) + (2. * cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p), 0.)));_T _acc{};
      { // subkernel 1
        const auto _interp22 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp30 = ZAcbc(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p2 * l1 * p + powr<2>(p)));
        const auto _interp31 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 + 2. * cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp32 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _interp33 = Zc(sqrt(powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _interp34 = RBdot(powr<2>(k), powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p));
        const auto _den12 = powr<-1>(_interp12 * _interp13 + _interp14 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den13 = powr<-2>(_interp13 * _interp32 + _interp33 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den14 = powr<-1>(_interp13 * _interp32 + _interp33 * (powr<2>(l1) + 2. * cosl1p2 * l1 * p + powr<2>(p)));
        const auto _den15 = powr<-1>(_interp13 * _interp26 + _interp29 * (powr<2>(l1) + 2. * (cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _cse1 = -50. * _interp13;
        const auto _cse2 = 50. * _interp28;_acc += fma(-0.25, _den12 * _den16 * _den2 * _den3 * _interp11 * _interp18 * _interp22 * (-1. + 2. * cosl1p1 * cosl1p2 + 2. * powr<2>(cosl1p2)) * (_interp3 * (-50. * _interp2 + 50. * _interp5) * powr<6>(k) + _interp1 * _interp2 * (1. + powr<6>(k)) + _interp3 * _interp4 * (1. + 1. * powr<6>(k))) * (2. * cosl1p1 * l1 + 4. * cosl1p2 * l1 - 3. * p) * p, fma(1., _den14 * _den17 * _den4 * _interp23 * ((_cse1 + _cse2) * _interp26 + _interp25 * _interp26 + _interp13 * _interp27) * _interp30 * _interp31 * (-0.5 + 1. * cosl1p1 * cosl1p2 + 1. * powr<2>(cosl1p2)) * p * (1. * cosl1p1 * l1 + 2. * cosl1p2 * l1 + 1.5 * p), fma(0.25, _den13 * _den15 * _den4 * _interp23 * _interp30 * _interp31 * ((_cse1 + _cse2 + _interp25) * _interp32 + _interp13 * _interp34) * (-1. + 2. * cosl1p1 * cosl1p2 + 2. * powr<2>(cosl1p2)) * p * (2. * cosl1p1 * l1 + 4. * cosl1p2 * l1 + 3. * p), 0.)));
      }
      { // subkernel 2
        const auto _interp7 = RB(powr<2>(k), powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _interp8 = ZA(sqrt(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _interp9 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) + cosl1p1 * l1 * p + powr<2>(p)));
        const auto _interp10 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (cosl1p1 - cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp15 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _interp16 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _interp17 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
        const auto _interp19 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) - 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _interp24 = ZAcbc(sqrt(0.6666666666666666 * powr<2>(l1) + 0.6666666666666666 * (2. * cosl1p1 + cosl1p2) * l1 * p + powr<2>(p)));
        const auto _den1 = powr<-1>(1. + powr<6>(k));
        const auto _den5 = powr<-1>(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _den6 = powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _den7 = powr<-1>(powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p));
        const auto _den8 = powr<-1>(_interp15 * _interp2 + _interp16 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den9 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den10 = powr<-1>(_interp2 * _interp7 + _interp8 * (powr<2>(l1) + 2. * cosl1p1 * l1 * p + powr<2>(p)));
        const auto _den11 = powr<-1>(_interp12 * _interp13 + _interp14 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));// clang-format off
        _acc += fma(3., _den10 * _den17 * _den4 * _den7 * _interp23 * _interp24 * (_interp25 * _interp26 + _interp13 * _interp27 + _interp26 * (-50. * _interp13 + 50. * _interp28)) * _interp9 * p * (0.6666666666666666 * cosl1p2 * powr<3>(l1) - 1.333333333333333 * powr<3>(cosl1p2) * powr<3>(l1) - powr<2>(l1) * p + 0.6666666666666666 * powr<2>(cosl1p2) * powr<2>(l1) * p - 0.3333333333333333 * cosl1p2 * l1 * powr<2>(p) - 0.5 * powr<3>(p) + powr<3>(cosl1p1) * l1 * p * (-0.6666666666666666 * cosl1p2 * l1 + 2.333333333333333 * p) + powr<2>(cosl1p1) * (-0.6666666666666666 * cosl1p2 * powr<3>(l1) + 2. * powr<2>(l1) * p - 2. * powr<2>(cosl1p2) * powr<2>(l1) * p + 3. * cosl1p2 * l1 * powr<2>(p) + 1. * powr<3>(p)) + cosl1p1 * ((0.3333333333333333 - 2. * powr<2>(cosl1p2)) * powr<3>(l1) + cosl1p2 * (2.333333333333333 - 1.333333333333333 * powr<2>(cosl1p2)) * powr<2>(l1) * p + (-1.166666666666667 + 0.6666666666666666 * powr<2>(cosl1p2)) * l1 * powr<2>(p) + 1. * cosl1p2 * powr<3>(p))), fma(0.5, _den16 * _den2 * _den3 * _den5 * _den8 * _interp17 * _interp18 * _interp19 * (_interp3 * (-50. * _interp2 + 50. * _interp5) * powr<6>(k) + _interp1 * _interp2 * (1. + powr<6>(k)) + _interp3 * _interp4 * (1. + 1. * powr<6>(k))) * p * (-4. * cosl1p2 * powr<3>(l1) + 8. * powr<3>(cosl1p2) * powr<3>(l1) - 6. * powr<2>(l1) * p + 4. * powr<2>(cosl1p2) * powr<2>(l1) * p + powr<3>(cosl1p1) * l1 * (-4. * cosl1p2 * l1 - 14. * p) * p + 2. * cosl1p2 * l1 * powr<2>(p) - 3. * powr<3>(p) + powr<2>(cosl1p1) * (4. * cosl1p2 * powr<3>(l1) + 12. * powr<2>(l1) * p - 12. * powr<2>(cosl1p2) * powr<2>(l1) * p - 18. * cosl1p2 * l1 * powr<2>(p) + 6. * powr<3>(p)) + cosl1p1 * ((-2. + 12. * powr<2>(cosl1p2)) * powr<3>(l1) + cosl1p2 * (14. - 8. * powr<2>(cosl1p2)) * powr<2>(l1) * p + (7. - 4. * powr<2>(cosl1p2)) * l1 * powr<2>(p) + 6. * cosl1p2 * powr<3>(p))), fma(-0.5, _den11 * _den3 * _den6 * _den9 * _interp10 * _interp11 * _interp9 * (_interp1 * _interp2 + _interp3 * (_interp4 + 50. * _den1 * (-_interp2 + _interp5) * powr<6>(k))) * p * (4. * cosl1p2 * powr<3>(l1) - 8. * powr<3>(cosl1p2) * powr<3>(l1) + 6. * powr<2>(l1) * p - 4. * powr<3>(cosl1p1) * cosl1p2 * powr<2>(l1) * p - 4. * powr<2>(cosl1p2) * powr<2>(l1) * p - 2. * cosl1p2 * l1 * powr<2>(p) + 3. * powr<3>(p) - 2. * powr<2>(cosl1p1) * l1 * (2. * cosl1p2 * powr<2>(l1) + l1 * p + 6. * powr<2>(cosl1p2) * l1 * p - 5. * cosl1p2 * powr<2>(p)) + cosl1p1 * ((2. - 12. * powr<2>(cosl1p2)) * powr<3>(l1) + (5. - 4. * powr<2>(cosl1p2)) * l1 * powr<2>(p) + 6. * cosl1p2 * powr<3>(p) + powr<2>(l1) * (6. * cosl1p2 * p - 8. * powr<3>(cosl1p2) * p))), 0.)));
        // clang-format on

      }
      return _acc;
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA)
    {
      using namespace DiFfRG;using namespace DiFfRG::compute;
      return 0.;
    }private: static KOKKOS_FORCEINLINE_FUNCTION auto RB(const auto& k2, const auto& p2)
    {
      return Regulator::RB(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto RF(const auto& k2, const auto& p2)
    {
      return Regulator::RF(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto RBdot(const auto& k2, const auto& p2)
    {
      return Regulator::RBdot(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto RFdot(const auto& k2, const auto& p2)
    {
      return Regulator::RFdot(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RB(const auto& k2, const auto& p2)
    {
      return Regulator::dq2RB(k2, p2);
    }

    static KOKKOS_FORCEINLINE_FUNCTION auto dq2RF(const auto& k2, const auto& p2)
    {
      return Regulator::dq2RF(k2, p2);
    }
  };
} using DiFfRG::ZAcbc_form_kernel;