#pragma once

#include "DiFfRG/common/utils.hh"
#include "DiFfRG/physics/integration.hh"
#include "./ZA/ZA.hh"
#include "./ZA3/ZA3.hh"
#include "./ZA4/ZA4.hh"
#include "./ZAAqbq1/ZAAqbq1.hh"
#include "./ZAAqbq2/ZAAqbq2.hh"
#include "./ZAcbc/ZAcbc.hh"
#include "./ZAqbq1/ZAqbq1.hh"
#include "./ZAqbq4/ZAqbq4.hh"
#include "./ZAqbq7/ZAqbq7.hh"
#include "./Zc/Zc.hh"
#include "./Zq/Zq.hh"

class QCDFlows
{
  public:
  QCDFlows(const DiFfRG::ConfigTree& config)
  ;

  void set_k(const double k)
  ;

  void set_T(const double T)
  ;

  void set_typical_E(const double E)
  ;

  void set_x_extent(const double x_extent)
  ;

  DiFfRG::QuadratureProvider quadrature_provider;

  ZA_integrator ZA;

  ZA3_integrator ZA3;

  ZA4_integrator ZA4;

  ZAAqbq1_integrator ZAAqbq1;

  ZAAqbq2_integrator ZAAqbq2;

  ZAcbc_integrator ZAcbc;

  ZAqbq1_integrator ZAqbq1;

  ZAqbq4_integrator ZAqbq4;

  ZAqbq7_integrator ZAqbq7;

  Zc_integrator Zc;

  Zq_integrator Zq;
};