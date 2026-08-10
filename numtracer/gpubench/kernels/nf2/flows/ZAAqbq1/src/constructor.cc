#include "../kernel.hh"

#include "../ZAAqbq1.hh"

ZAAqbq1_integrator::ZAAqbq1_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config) : integrator(quadrature_provider, config), quadrature_provider(quadrature_provider)
{
}