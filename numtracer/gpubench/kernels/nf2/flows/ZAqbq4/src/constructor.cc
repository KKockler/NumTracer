#include "../kernel.hh"

#include "../ZAqbq4.hh"

ZAqbq4_integrator::ZAqbq4_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config) : integrator(quadrature_provider, config), quadrature_provider(quadrature_provider)
{
}