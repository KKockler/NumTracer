#pragma once

#include "shim.hpp"

namespace DiFfRG {
  template<typename _Regulator>
  class EtaPiL_form_kernel
  {
    public:
    using Regulator = _Regulator;

    // clang-format off
static KOKKOS_FORCEINLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
{
using namespace DiFfRG;
    // clang-format on
    using namespace DiFfRG::compute;const auto _interp1 = RFdot(powr<2>(k), powr<2>(l1));
    const auto _interp2 = Zq(k);
    const auto _interp3 = RF(powr<2>(k), powr<2>(l1));
    const auto _interp4 = dtZq(k);
    const auto _interp5 = Zq(1.02 * k);
    const auto _interp6 = hSigL(l1);
    const auto _interp7 = Zq(l1);
    const auto _interp8 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
    const auto _interp9 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
    const auto _interp10 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
    const auto _interp11 = hPiL(0.000816496580927726 * sqrt(1. - 1000. * (cos1 - 1000. * l1) * l1));
    const auto _interp12 = hPiL(0.816496580927726 * l1);
    const auto _den1 = powr<-3>(powr<2>(_interp2 * _interp3 + _interp7 * l1) + 0.5 * powr<2>(_interp6) * rhoL);
    const auto _den2 = powr<-1>(powr<2>(_interp2) * powr<2>(_interp9) + 0.002 * _interp10 * _interp2 * _interp9 * sqrt(1. - 2000. * (cos1 - 500. * l1) * l1) + powr<2>(_interp10) * (1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)) + 0.5 * powr<2>(_interp8) * rhoL);
    const auto _cse1 = -0.002 * cos1 * l1;
    const auto _cse2 = 1.e-6 * fma(1.e6, _cse1, fma(1.e6, powr<2>(l1), 1.));
    const auto _cse3 = sqrt(_cse2);
    const auto _cse4 = _interp7 * l1;
    const auto _cse5 = _interp2 * _interp3;
    const auto _cse6 = _cse4 + _cse5;
    const auto _cse7 = sqrt(powr<-1>(_cse2));
    const auto _cse8 = _cse3 * _interp10;
    const auto _cse9 = _interp2 * _interp9;
    const auto _cse10 = _cse8 + _cse9;
    return -9.e6 * _den1 * _den2 * fma(-1., _interp1 * _interp2, fma(-1., _interp3 * (-50. * _interp2 + _interp4 + 50. * _interp5), 0.)) * fma(-0.001, _cse10 * _cse7 * powr<2>(_interp11) * powr<4>(_interp6) * cos1 * powr<2>(rhoL), fma(2., _cse6 * powr<2>(_interp12) * (powr<2>(_cse6) + 0.5 * powr<2>(_interp6) * rhoL) * (2. * _cse2 * powr<2>(_interp10) + 4. * _cse3 * _interp10 * _interp2 * _interp9 + 2. * powr<2>(_interp2) * powr<2>(_interp9) + powr<2>(_interp8) * rhoL), fma(-1., _cse7 * powr<2>(_interp11) * (4. * _cse10 * powr<4>(_cse6) * l1 + powr<3>(_interp6) * (2. * _cse3 * _cse6 * _interp8 - _cse10 * _interp6 * l1) * powr<2>(rhoL) + 4. * powr<3>(_cse6) * (-0.001 * _cse10 * _cse6 * cos1 + _cse3 * _interp6 * _interp8 * rhoL)), 0.)));
  }

  // clang-format off
static KOKKOS_FORCEINLINE_FUNCTION auto constant(const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
{
using namespace DiFfRG;
  // clang-format on
  using namespace DiFfRG::compute;
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
} using DiFfRG::EtaPiL_form_kernel;