#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Colsum_num_kernels.hh"

namespace numtracer_kernels
{
  class Colsum_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& cos2, const double& p)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::colsum_num::nenv) > 0 ? (numtracer_kernels::colsum_num::nenv) : 1];
      numtracer_kernels::colsum_num::fill(fenv, l1, cos1, cos2, p);  return fma(-2., numtracer_kernels::colsum_num::tr1(fenv), numtracer_kernels::colsum_num::tr0(fenv));
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Colsum_num_kernel;