#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Epsstraddle_num_kernels.hh"

namespace numtracer_kernels
{
  class Epsstraddle_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::epsstraddle_num::nenv) > 0 ? (numtracer_kernels::epsstraddle_num::nenv) : 1];
      numtracer_kernels::epsstraddle_num::fill(fenv, l1, cos1, cos2, p);  return fma(2., numtracer_kernels::epsstraddle_num::tr0(fenv), fma(-2., numtracer_kernels::epsstraddle_num::tr1(fenv), fma(3., numtracer_kernels::epsstraddle_num::tr2(fenv), fma(-3., numtracer_kernels::epsstraddle_num::tr3(fenv), 0.))));
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Epsstraddle_num_kernel;