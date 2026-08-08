#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Colpow_num_kernels.hh"

namespace numtracer_kernels
{
  class Colpow_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::colpow_num::nenv) > 0 ? (numtracer_kernels::colpow_num::nenv) : 1];
      numtracer_kernels::colpow_num::fill(fenv, l1, cos1, cos2, p);
      return numtracer_kernels::colpow_num::tr0(fenv);
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Colpow_num_kernel;