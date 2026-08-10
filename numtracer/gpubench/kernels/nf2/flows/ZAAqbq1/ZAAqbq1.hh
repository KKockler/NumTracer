#pragma once

#include "DiFfRG/physics/integration.hh"
#include "DiFfRG/physics/physics.hh"
#include "DiFfRG/physics/interpolation.hh"

namespace DiFfRG { template<typename> class ZAAqbq1_kernel;

  class ZAAqbq1_integrator
  {
    public:
    ZAAqbq1_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config)
    ;


    using Regulator = DiFfRG::PolynomialExpRegulator<>;

    Integrator_p2_4D_3ang<4, double, ZAAqbq1_kernel<Regulator>, DiFfRG::GPU_exec> integrator;

    // clang-format off
    DiFfRG::GPU_exec map(double* dest, const LogarithmicCoordinates1D<double>& coordinates, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2)
    // clang-format on
    ;

    template<typename IT, typename C, typename ...T>
    DiFfRG::GPU_exec map(IT* dest, const C& coordinates, const device::tuple<T...>& args)
    {
      return device::apply([&](const auto&...t){return map(dest, coordinates, t...);}, args);
    }

    // clang-format off
    void get(double& dest, const double& p, const double& k, const double& etaQ, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA3, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAcbc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq4, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAqbq7, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& Zc, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& dtZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZA, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& zq, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq1, const SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>& ZAAqbq2)
    // clang-format on
    ;

    template<typename IT, typename ...T>
    void get(IT& dest, const double& p, const device::tuple<T...>& args)
    {
      device::apply([&](const auto&...t){get(dest, p, t...);}, args);
    }
    private:
    DiFfRG::QuadratureProvider& quadrature_provider;
  };
}
using DiFfRG::ZAAqbq1_integrator;