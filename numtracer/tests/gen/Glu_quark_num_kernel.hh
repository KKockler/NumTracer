#pragma once

#include "nt_regulators.hpp"
#include "numtracer/codegen/runtime.hpp"
#include "numtracer/sun/sun_data.hpp"
#include "Glu_quark_num_kernels.hh"

namespace numtracer_kernels
{
  class Glu_quark_num_kernel
  {
    public:
    static inline auto kernel(const double& l1, const double& cos1, const double& p, const double& k)
    {
      using namespace numtracer;
      using namespace numtracer::compute;
      double fenv[(numtracer_kernels::glu_quark_num::nenv) > 0 ? (numtracer_kernels::glu_quark_num::nenv) : 1];
      numtracer_kernels::glu_quark_num::fill(fenv, l1, cos1, p);
      const auto _interp1 = RF(powr<2>(k), powr<2>(l1));
      const auto _interp2 = RF(powr<2>(k), powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      const auto _interp3 = RFdot(powr<2>(k), powr<2>(l1));
      const auto _den1 = powr<-2>(_interp1 + powr<2>(l1));
      const auto _den2 = powr<-1>(_interp2 + powr<2>(l1) - 2. * cos1 * l1 * p + powr<2>(p));
      return -0.1666666666666667 * _den1 * _den2 * _interp3 * numtracer_kernels::glu_quark_num::tr0(fenv);
    }

    static inline auto constant(const double& p, const double& k)
    {
      return 0.;
    }
  };
}
using numtracer_kernels::Glu_quark_num_kernel;