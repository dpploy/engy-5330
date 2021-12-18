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

#include "NormalFluxBC.h"
#include "Function.h"

registerMooseObject("FuelPlateApp", NormalFluxBC);

defineLegacyParams(NormalFluxBC);

InputParameters
NormalFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Robin type normal flux boundary condition.");
  params.addParam<Real>("transferCoeff", 0.0, "Transfer coefficient");
  return params;
}

NormalFluxBC::NormalFluxBC(const InputParameters & parameters): 
  IntegratedBC(parameters),
  _transferCoeff(getParam<Real>("transferCoeff"))
{
}

Real
NormalFluxBC::computeQpResidual()
{
 // Residual
 return _transferCoeff * _u[_qp] * _test[_i][_qp];
}

Real
NormalFluxBC::computeQpJacobian()
{
 // Jacobian
 return _transferCoeff * _phi[_j][_qp] * _test[_i][_qp];
}
