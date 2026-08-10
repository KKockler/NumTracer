#include "../kernel.hh"

#include "../ZAAqbq2.hh"

ZAAqbq2_integrator::ZAAqbq2_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config) : integrator(quadrature_provider, config), quadrature_provider(quadrature_provider)
{
}