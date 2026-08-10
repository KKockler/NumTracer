#include "../kernel.hh"

#include "../ZAqbq1.hh"

ZAqbq1_integrator::ZAqbq1_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config) : integrator(quadrature_provider, config), quadrature_provider(quadrature_provider)
{
}