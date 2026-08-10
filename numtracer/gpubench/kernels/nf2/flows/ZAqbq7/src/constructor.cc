#include "../kernel.hh"

#include "../ZAqbq7.hh"

ZAqbq7_integrator::ZAqbq7_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config) : integrator(quadrature_provider, config), quadrature_provider(quadrature_provider)
{
}