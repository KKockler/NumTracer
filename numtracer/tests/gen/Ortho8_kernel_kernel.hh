#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Ortho8_kernel_kernels.hh"

namespace numtracer_kernels
{
  template<typename REG>
  class Ortho8_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::ortho8::nenv) > 0 ? (numtracer_kernels::ortho8::nenv) : 1];
      numtracer_kernels::ortho8::fill(fenv, l1, cos1, cos2, p);  return 2. * numtracer_kernels::ortho8::tr0(fenv) * Global_TBGetInnerProduct(AqbqDirect8)(Global_TBGetProjector, 1., Global_TBGetBasisElement, 1.);
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
using numtracer_kernels::Ortho8_kernel;