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

    static KOKKOS_INLINE_FUNCTION auto kernel(const double& l1, const double& cos1, const double& cos2, const double& phi, const double& p, const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    {
      using namespace DiFfRG;
      using namespace DiFfRG::compute;
      using namespace numtracer;
      double fenv[(DiFfRG::za4_ym::nenv) > 0 ? (DiFfRG::za4_ym::nenv) : 1];
      DiFfRG::za4_ym::fill(fenv, l1, cos1, cos2, phi, p);
      const double cosl1p1 = sqrt(1. - powr<2>(cos1)) * cos2;
      const double cosl1p2 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 - 2. * sqrt(2. - 2. * powr<2>(cos2)) * cos(phi));
      const double cosl1p3 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) - 1.732050807568877 * sin(phi)));
      const double cosl1p4 = -0.3333333333333333 * sqrt(1. - powr<2>(cos1)) * (cos2 + sqrt(2. - 2. * powr<2>(cos2)) * (cos(phi) + 1.732050807568877 * sin(phi)));
      const auto _interp1 = RB(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RBdot(powr<2>(k), powr<2>(l1));
      const auto _ix4 = ZA.index(l1);
      const auto _interp3 = ZA.at(_ix4);
      const auto _interp4 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p));
      const auto _ix7 = ZA.index(sqrt(powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp5 = ZA.at(_ix7);
      const auto _interp6 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p));
      const auto _ix5 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp7 = ZA.at(_ix5);
      const auto _interp8 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _ix6 = ZA.index(sqrt(powr<2>(l1) - 2. * (cosl1p1 + cosl1p2) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp9 = ZA.at(_ix6);
      const auto _ix3 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p1 * l1 * p + powr<2>(p)));
      const auto _interp10 = ZA3.at(_ix3);
      const auto _ix2 = ZA3.index(0.816496580927726 * sqrt(powr<2>(l1) - (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _interp11 = ZA3.at(_ix2);
      const auto _ix1 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * cosl1p1 + cosl1p2) * l1 * p + 10. * powr<2>(p)));
      const auto _interp12 = ZA3.at(_ix1);
      const auto _ix0 = ZA3.index(0.3333333333333333 * sqrt(6. * powr<2>(l1) - 6. * (2. * (cosl1p1 + cosl1p2) + cosl1p3) * l1 * p + 10. * powr<2>(p)));
      const auto _interp13 = ZA3.at(_ix0);
      const auto _interp14 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p1 + cosl1p2) * l1 * p + 5. * powr<2>(p)));
      const auto _interp15 = RB(powr<2>(k), powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p));
      const auto _interp16 = ZA(sqrt(powr<2>(l1) - 2. * (cosl1p2 + cosl1p3) * l1 * p + 1.333333333333333 * powr<2>(p)));
      const auto _interp17 = ZA4(0.408248290463863 * sqrt(3. * powr<2>(l1) - 3. * (cosl1p2 + cosl1p3) * l1 * p + 5. * powr<2>(p)));
      const auto _interp18 = RB(powr<2>(k), powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p));
      const auto _interp19 = ZA(sqrt(powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp20 = ZA3(0.816496580927726 * sqrt(powr<2>(l1) - cosl1p2 * l1 * p + powr<2>(p)));
      const auto _interp21 = ZA4(0.7071067811865475 * sqrt(powr<2>(l1) - (cosl1p1 + 2. * cosl1p2 + cosl1p3) * l1 * p + 2. * powr<2>(p)));
      const auto _interp22 = ZAcbc.at(_ix3);
      const auto _interp23 = ZAcbc.at(_ix2);
      const auto _interp24 = ZAcbc.at(_ix1);
      const auto _interp25 = ZAcbc.at(_ix0);
      const auto _interp26 = Zc.at(_ix4);
      const auto _interp27 = Zc.at(_ix7);
      const auto _interp28 = Zc.at(_ix5);
      const auto _interp29 = Zc.at(_ix6);
      const auto _den1 = powr<-1>(1. + powr<6>(k));
      const auto _den2 = powr<-2>(_interp3 * powr<2>(l1) + _interp1 * nthk0);
      const auto _den3 = powr<-2>(_interp26 * powr<2>(l1) + _interp1 * nthk3);
      const auto _den4 = powr<-1>(_interp4 * nthk0 + _interp5 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den5 = powr<-1>(_interp18 * nthk0 + _interp19 * (powr<2>(l1) - 2. * cosl1p2 * l1 * p + powr<2>(p)));
      const auto _den6 = -powr<-1>(_interp4 * nthk3 + _interp27 * (powr<2>(l1) - 2. * cosl1p1 * l1 * p + powr<2>(p)));
      const auto _den7 = 3. * powr<-1>(3. * _interp8 * nthk0 + _interp9 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den8 = 3. * powr<-1>(3. * _interp15 * nthk0 + _interp16 * (3. * powr<2>(l1) - 6. * (cosl1p2 + cosl1p3) * l1 * p + 4. * powr<2>(p)));
      const auto _den9 = powr<-1>(_interp6 * nthk0 + _interp7 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _den10 = -3. * powr<-1>(3. * _interp8 * nthk3 + _interp29 * (3. * powr<2>(l1) - 6. * (cosl1p1 + cosl1p2) * l1 * p + 4. * powr<2>(p)));
      const auto _den11 = -powr<-1>(_interp6 * nthk3 + _interp28 * (powr<2>(l1) - 2. * (cosl1p1 + cosl1p2 + cosl1p3) * l1 * p + powr<2>(p)));
      const auto _cse1 = _interp2 * nthk0;
      const auto _cse2 = -nthk0;
      const auto _cse3 = _cse2 + nthk2;
      const auto _cse4 = 50. * _cse3 * _den1 * powr<6>(k);
      const auto _cse5 = _cse4 + nthk1;
      const auto _cse6 = _cse5 * _interp1;
      const auto _cse7 = _cse1 + _cse6;
      return fma(0.00510204081632653, _cse7 * _den2 * _den4 * _den7 * _den9 * DiFfRG::za4_ym::tr0(fenv) * _interp10 * _interp11 * _interp12 * _interp13, fma(0.00510204081632653, _cse7 * _den2 * _den7 * _den9 * DiFfRG::za4_ym::tr1(fenv) * _interp11 * _interp13 * _interp14, fma(0.001275510204081632, _cse7 * _den2 * _den8 * DiFfRG::za4_ym::tr3(fenv) * powr<2>(_interp17), fma(0.002551020408163265, _cse7 * _den2 * _den5 * _den9 * DiFfRG::za4_ym::tr2(fenv) * _interp11 * _interp20 * _interp21, fma(0.01020408163265306, _den10 * _den11 * _den3 * _den6 * DiFfRG::za4_ym::tr4(fenv) * _interp22 * _interp23 * _interp24 * _interp25 * (_interp2 * nthk3 + _interp1 * (nthk4 + 50. * (-nthk3 + nthk5))), 0.)))));
    }

    static KOKKOS_INLINE_FUNCTION auto constant(const double& p, const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA, const double& nthk0, const double& nthk1, const double& nthk2, const double& nthk3, const double& nthk4, const double& nthk5)
    {
      return 0.;
    }

    static device::array<double, 6> ntHoisted(const double& k, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, FocusedLogCoordinates1D<double>, GPU_memory>& ZA)
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