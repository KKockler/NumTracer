#pragma once

#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "EMVert_diag_kernels.hh"

namespace numtracer_kernels
{
  class EMVert_diag_kernel
  {
    public:
    static inline auto kernel(const double& p01, const double& p02, const double& p03, const double& p, const double& f0, const double& l1, const double& cos1, const double& phi)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::emvert_diag::nenv) > 0 ? (numtracer_kernels::emvert_diag::nenv) : 1];
      numtracer_kernels::emvert_diag::fill(fenv, p01, p02, p03, p, f0, l1, cos1, phi);
      return numtracer_kernels::emvert_diag::tr0(fenv);
    }

    static inline auto constant(const double& p)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::EMVert_diag_kernel;