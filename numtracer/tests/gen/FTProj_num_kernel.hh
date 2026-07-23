#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "FTProj_num_kernels.hh"

namespace numtracer_kernels
{
  class FTProj_num_kernel
  {
    public:
    static inline auto kernel(const double& p0, const double& p, const double& l0, const double& l1, const double& cos1)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::ftproj_num::nenv) > 0 ? (numtracer_kernels::ftproj_num::nenv) : 1];
      numtracer_kernels::ftproj_num::fill(fenv, p0, p, l0, l1, cos1);  return fma(numtracer_kernels::ftproj_num::tr1(fenv), cos1 * l0 * l1 * p, numtracer_kernels::ftproj_num::tr0(fenv));
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::FTProj_num_kernel;