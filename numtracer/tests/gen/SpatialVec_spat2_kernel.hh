#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "SpatialVec_spat2_kernels.hh"

namespace numtracer_kernels
{
  class SpatialVec_spat2_kernel
  {
    public:
    static inline auto kernel(const double& p0, const double& p, const double& l0, const double& l1, const double& cos1)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::spatialvec_spat2::nenv) > 0 ? (numtracer_kernels::spatialvec_spat2::nenv) : 1];
      numtracer_kernels::spatialvec_spat2::fill(fenv, p0, p, l0, l1, cos1);
      return -numtracer_kernels::spatialvec_spat2::tr0(fenv);
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::SpatialVec_spat2_kernel;