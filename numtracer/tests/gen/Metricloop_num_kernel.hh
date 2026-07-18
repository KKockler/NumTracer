#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Metricloop_num_kernels.hh"

namespace numtracer_kernels
{
  template<typename REG>
  class Metricloop_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::metricloop_num::nenv) > 0 ? (numtracer_kernels::metricloop_num::nenv) : 1];
      numtracer_kernels::metricloop_num::fill(fenv, l1, cos1, cos2, p);  return numtracer_kernels::metricloop_num::tr0(fenv);
    }

    static inline auto constant(const double& p)
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
using numtracer_kernels::Metricloop_num_kernel;