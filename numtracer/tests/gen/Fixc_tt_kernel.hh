#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Fixc_tt_kernels.hh"

namespace numtracer_kernels
{
  template<typename REG>
  class Fixc_tt_kernel
  {
    public:
    static inline auto kernel(const double& p0, const double& p, const double& l0, const double& l1, const double& cos1)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::fixc_tt::nenv) > 0 ? (numtracer_kernels::fixc_tt::nenv) : 1];
      numtracer_kernels::fixc_tt::fill(fenv, p0, p, l0, l1, cos1);  return numtracer_kernels::fixc_tt::tr0(fenv);
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
using numtracer_kernels::Fixc_tt_kernel;