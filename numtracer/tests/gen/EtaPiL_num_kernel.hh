#pragma once

#include "shim.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "EtaPiL_num_kernels.hh"

namespace DiFfRG
{
  class EtaPiL_num_kernel
  {
    public:
    // clang-format off
static inline auto kernel(const double& l1, const double& cos1, const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
{
using namespace DiFfRG;
    // clang-format on

    using namespace DiFfRG::compute;
    using namespace numtracer;
    double fenv[(DiFfRG::etapil_num::nenv) > 0 ? (DiFfRG::etapil_num::nenv) : 1];
    const double dr_0 = 0.7071067811865475;
    const double dr_1 = sqrt(rhoL);
    const double dr_2 = hSigL(l1);
    const double dr_3 = -powr<-1>(l1) * RF(powr<2>(k), powr<2>(l1)) * Zq(k) - Zq(l1);
    DiFfRG::etapil_num::fill(fenv, l1, cos1, dr_0, dr_1, dr_2, dr_3);
    const auto _interp1 = hPiL(0.816496580927726 * l1);
    const auto _interp2 = ntRe(DiFfRG::etapil_num::tr0(fenv));
    const auto _interp3 = RFdot(powr<2>(k), powr<2>(l1));
    const auto _interp4 = Zq(k);
    const auto _interp5 = RF(powr<2>(k), powr<2>(l1));
    const auto _interp6 = dtZq(k);
    const auto _interp7 = Zq(1.02 * k);
    const auto _interp8 = hSigL(l1);
    const auto _interp9 = Zq(l1);
    const auto _interp10 = ntIm(DiFfRG::etapil_num::tr1(fenv));
    const auto _interp11 = hPiL(0.000816496580927726 * sqrt(1. - 1000. * (cos1 - 1000. * l1) * l1));
    const auto _interp12 = hSigL(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
    const auto _interp13 = RF(powr<2>(k), 1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1));
    const auto _interp14 = Zq(sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1)));
    const auto _interp15 = ntRe(DiFfRG::etapil_num::tr2(fenv));
    const auto _interp16 = ntRe(DiFfRG::etapil_num::tr3(fenv));
    const auto _den1 = powr<-3>(powr<2>(_interp4 * _interp5 + _interp9 * l1) + 0.5 * powr<2>(_interp8) * rhoL);
    const auto _den2 = powr<-1>(powr<2>(_interp13 * _interp4 + _interp14 * sqrt(1.e-6 - 0.002 * cos1 * l1 + powr<2>(l1))) + 0.5 * powr<2>(_interp12) * rhoL);
    const auto _cse1 = powr<-1>(l1);
    const auto _cse2 = -_interp9;
    const auto _cse3 = -_cse1 * _interp4 * _interp5;
    const auto _cse4 = _cse2 + _cse3;
    const auto _cse5 = -_cse1 * _interp3 * _interp4;
    const auto _cse6 = -_interp4;
    const auto _cse7 = _cse6 + _interp7;
    const auto _cse8 = 50. * _cse7;
    const auto _cse9 = _cse8 + _interp6;
    const auto _cse10 = -_cse1 * _cse9 * _interp5;
    const auto _cse11 = _cse10 + _cse5;
    const auto _cse12 = -_interp14;
    const auto _cse13 = -0.002 * cos1 * l1;
    const auto _cse14 = 1.e-6 * fma(1.e6, _cse13, fma(1.e6, powr<2>(l1), 1.));
    const auto _cse15 = sqrt(powr<-1>(_cse14));
    const auto _cse16 = -_cse15 * _interp13 * _interp4;
    const auto _cse17 = _cse12 + _cse16;
    const auto _cse18 = sqrt(rhoL);
    return 500000. * fma(4., _cse11 * _cse4 * _den1 * powr<2>(_interp1) * _interp2, fma(-2., 1.414213562373095 * _cse11 * _cse14 * powr<2>(_cse17) * _cse18 * _den1 * _den2 * powr<2>(_interp1) * _interp10 * _interp8, fma(2., 1.414213562373095 * _cse11 * _cse18 * powr<2>(_cse4) * _den1 * _den2 * _interp10 * powr<2>(_interp11) * _interp12 * powr<2>(l1), fma(-4., _cse11 * _cse17 * powr<2>(_cse4) * _den1 * _den2 * powr<2>(_interp11) * _interp16 * powr<2>(l1), fma(-1., 1.414213562373095 * _cse11 * _cse18 * _den1 * _den2 * powr<2>(_interp1) * _interp10 * powr<2>(_interp12) * _interp8 * rhoL, fma(-2., _cse11 * _cse17 * _den1 * _den2 * powr<2>(_interp11) * _interp15 * powr<2>(_interp8) * rhoL, fma(1.414213562373095, _cse11 * _den1 * _den2 * _interp10 * powr<2>(_interp11) * _interp12 * powr<2>(_interp8) * sqrt(powr<3>(rhoL)), 0.)))))));
  }

  // clang-format off
static inline auto constant(const double& k, const double& etaPiL, const double& etaSigL, const double& d1V, const double& d2V, const double& rhoL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& hSigL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zPiL, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zSigL)
{
return 0.;
  // clang-format on

}
private:
static inline double ntRe(double x) { return x; }
template <class T> static inline double ntRe(const T &z) { return z.real(); }
static inline double ntIm(double) { return 0.0; }
template <class T> static inline double ntIm(const T &z) { return z.imag(); }
};
}
using DiFfRG::EtaPiL_num_kernel;