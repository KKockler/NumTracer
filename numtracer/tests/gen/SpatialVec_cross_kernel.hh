#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "SpatialVec_cross_kernels.hh"

namespace numtracer_kernels
{
  class SpatialVec_cross_kernel
  {
    public:
    static inline auto kernel(const double& p0, const double& p, const double& l0, const double& l1, const double& cos1)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::spatialvec_cross::nenv) > 0 ? (numtracer_kernels::spatialvec_cross::nenv) : 1];
      numtracer_kernels::spatialvec_cross::fill(fenv, p0, p, l0, l1, cos1);
      return fma(numtracer_kernels::spatialvec_cross::tr1(fenv), numtracer_kernels::spatialvec_cross::tr2(fenv) * numtracer_kernels::spatialvec_cross::tr3(fenv), fma(-4., numtracer_kernels::spatialvec_cross::tr0(fenv) * cos1 * l1 * p, 0.));
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::SpatialVec_cross_kernel;