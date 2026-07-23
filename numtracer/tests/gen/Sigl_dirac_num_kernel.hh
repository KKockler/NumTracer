#pragma once

#include "nt_regulators.hpp"
#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Sigl_dirac_num_kernels.hh"

namespace numtracer_kernels
{
  class Sigl_dirac_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::sigl_dirac_num::nenv) > 0 ? (numtracer_kernels::sigl_dirac_num::nenv) : 1];
      numtracer_kernels::sigl_dirac_num::fill(fenv, l1, cos1, p);
      const auto _interp1 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RF(powr<2>(k), powr<2>(p));
      const auto _den1 = powr<-1>(_interp2 + powr<2>(p));
      const auto _den2 = powr<-1>(_interp1 + powr<2>(l1));
      return -0.25 * _den1 * _den2 * numtracer_kernels::sigl_dirac_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const double& k)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Sigl_dirac_num_kernel;