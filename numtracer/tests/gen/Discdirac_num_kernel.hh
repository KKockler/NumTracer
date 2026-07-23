#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Discdirac_num_kernels.hh"

namespace numtracer_kernels
{
  class Discdirac_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::discdirac_num::nenv) > 0 ? (numtracer_kernels::discdirac_num::nenv) : 1];
      numtracer_kernels::discdirac_num::fill(fenv, l1, cos1, cos2, p);  return fma(numtracer_kernels::discdirac_num::tr0(fenv), numtracer_kernels::discdirac_num::tr3(fenv) * numtracer_kernels::discdirac_num::tr4(fenv), fma(numtracer_kernels::discdirac_num::tr1(fenv), numtracer_kernels::discdirac_num::tr5(fenv) * numtracer_kernels::discdirac_num::tr6(fenv), fma(numtracer_kernels::discdirac_num::tr2(fenv), numtracer_kernels::discdirac_num::tr7(fenv) * numtracer_kernels::discdirac_num::tr8(fenv), 0.)));
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Discdirac_num_kernel;