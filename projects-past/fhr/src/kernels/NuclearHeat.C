//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Engy-5310: Computational Continuum Transport Phenomena
//* UMass Lowell, Nuclear Chemical Engineering
//* https://github.com/dpploy/engy-5310

#include "NuclearHeat.h"

registerMooseObject("FHRApp", NuclearHeat);

template<>
InputParameters validParams<NuclearHeat>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Nuclear heat source term");
  params.addParam<Real>("sourceS", 1.0, "Nuclear heating power density");
  return params;
}

NuclearHeat::NuclearHeat(const InputParameters & parameters): 
  Kernel(parameters),
  sourceS(getParam<Real>("sourceS"))
{
}

Real
NuclearHeat::computeQpResidual()
{
  // Residual
  return - sourceS * _test[_i][_qp];
}

Real
NuclearHeat::computeQpJacobian()
{
  // Jacobian diagonal
  return 0.0;
}
