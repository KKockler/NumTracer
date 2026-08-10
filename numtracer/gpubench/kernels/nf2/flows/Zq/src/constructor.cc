#include "../kernel.hh"

#include "../Zq.hh"

Zq_integrator::Zq_integrator(DiFfRG::QuadratureProvider& quadrature_provider, const DiFfRG::ConfigTree& config) : integrator(quadrature_provider, config), quadrature_provider(quadrature_provider)
{
}