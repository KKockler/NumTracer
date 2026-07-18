#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Ortho147_kernels.hh"

namespace numtracer_kernels
{
  template<typename REG>
  class Ortho147_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p1m, const double& p2m, const double& cosP)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::ortho147::nenv) > 0 ? (numtracer_kernels::ortho147::nenv) : 1];
      numtracer_kernels::ortho147::fill(fenv, l1, cos1, cos2, p1m, p2m, cosP);
      const auto _den1 = powr<-2>(p1m);
      const auto _den2 = powr<-1>(-1. + powr<2>(cosP)) * powr<-2>(p1m) * powr<-2>(p2m);
      return 0.005208333333333332 * fma(8., numtracer_kernels::ortho147::tr0(fenv), fma(10., numtracer_kernels::ortho147::tr1(fenv), fma(-44., _den1 * numtracer_kernels::ortho147::tr2(fenv), fma(-102., _den2 * numtracer_kernels::ortho147::tr3(fenv), fma(-69., _den2 * numtracer_kernels::ortho147::tr4(fenv), 0.)))));
    }

    static inline auto constant()
    {
      return 0.;
    }
    private:
    static inline auto RB(const auto &k2, const auto &p2) { return REG::RB(k2, p2); }
    static inline auto RF(const auto &k2, const auto &p2) { return REG::RF(k2, p2); }
    static inline auto RBdot(const auto &k2, const auto &p2) { return REG::RBdot(k2, p2); }
    static inline auto RFdot(const auto &k2, const auto &p2) { return REG::RFdot(k2, p2); }
    static inline auto dq2RB(const auto &k2, const auto &p2) { return REG::dq2RB(k2, p2); }
    static inline auto dq2RF(const auto &k2, const auto &p2) { return REG::dq2RF(k2, p2); }
  };
}
using numtracer_kernels::Ortho147_kernel;