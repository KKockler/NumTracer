#include "./flows.hh"

QCDFlows::QCDFlows(const DiFfRG::ConfigTree& config) : quadrature_provider(config), ZA(quadrature_provider, config), ZA3(quadrature_provider, config), ZA4(quadrature_provider, config), ZAAqbq1(quadrature_provider, config), ZAAqbq2(quadrature_provider, config), ZAcbc(quadrature_provider, config), ZAqbq1(quadrature_provider, config), ZAqbq4(quadrature_provider, config), ZAqbq7(quadrature_provider, config), Zc(quadrature_provider, config), Zq(quadrature_provider, config)
{
}
void QCDFlows::set_k(const double k)
{
  DiFfRG::all_set_k(ZA, k);DiFfRG::all_set_k(ZA3, k);DiFfRG::all_set_k(ZA4, k);DiFfRG::all_set_k(ZAAqbq1, k);DiFfRG::all_set_k(ZAAqbq2, k);DiFfRG::all_set_k(ZAcbc, k);DiFfRG::all_set_k(ZAqbq1, k);DiFfRG::all_set_k(ZAqbq4, k);DiFfRG::all_set_k(ZAqbq7, k);DiFfRG::all_set_k(Zc, k);DiFfRG::all_set_k(Zq, k);
}
void QCDFlows::set_T(const double T)
{
  DiFfRG::all_set_T(ZA, T);DiFfRG::all_set_T(ZA3, T);DiFfRG::all_set_T(ZA4, T);DiFfRG::all_set_T(ZAAqbq1, T);DiFfRG::all_set_T(ZAAqbq2, T);DiFfRG::all_set_T(ZAcbc, T);DiFfRG::all_set_T(ZAqbq1, T);DiFfRG::all_set_T(ZAqbq4, T);DiFfRG::all_set_T(ZAqbq7, T);DiFfRG::all_set_T(Zc, T);DiFfRG::all_set_T(Zq, T);
}
void QCDFlows::set_typical_E(const double E)
{
  DiFfRG::all_set_typical_E(ZA, E);DiFfRG::all_set_typical_E(ZA3, E);DiFfRG::all_set_typical_E(ZA4, E);DiFfRG::all_set_typical_E(ZAAqbq1, E);DiFfRG::all_set_typical_E(ZAAqbq2, E);DiFfRG::all_set_typical_E(ZAcbc, E);DiFfRG::all_set_typical_E(ZAqbq1, E);DiFfRG::all_set_typical_E(ZAqbq4, E);DiFfRG::all_set_typical_E(ZAqbq7, E);DiFfRG::all_set_typical_E(Zc, E);DiFfRG::all_set_typical_E(Zq, E);
}
void QCDFlows::set_x_extent(const double x_extent)
{
  DiFfRG::all_set_x_extent(ZA, x_extent);DiFfRG::all_set_x_extent(ZA3, x_extent);DiFfRG::all_set_x_extent(ZA4, x_extent);DiFfRG::all_set_x_extent(ZAAqbq1, x_extent);DiFfRG::all_set_x_extent(ZAAqbq2, x_extent);DiFfRG::all_set_x_extent(ZAcbc, x_extent);DiFfRG::all_set_x_extent(ZAqbq1, x_extent);DiFfRG::all_set_x_extent(ZAqbq4, x_extent);DiFfRG::all_set_x_extent(ZAqbq7, x_extent);DiFfRG::all_set_x_extent(Zc, x_extent);DiFfRG::all_set_x_extent(Zq, x_extent);
}