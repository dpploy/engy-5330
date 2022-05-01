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

#include "NormalHeatFluxBC.h"
#include "Function.h"

registerMooseObject("HeatedPlateApp", NormalHeatFluxBC);

defineLegacyParams(NormalHeatFluxBC);

InputParameters
NormalHeatFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addClassDescription("Robin type heat normal flux boundary condition.");
  params.addRequiredParam<FunctionName>("refTempFunc", "Profile function");
  params.addParam<Real>("transferCoeff", 0.0, "Heat transfer coefficient");
  return params;
}

NormalHeatFluxBC::NormalHeatFluxBC(const InputParameters & parameters): 
  IntegratedBC(parameters),
  _refTempFunc(getFunction("refTempFunc")),
  _transferCoeff(getParam<Real>("transferCoeff"))
{
}

Real
NormalHeatFluxBC::computeQpResidual()
{
 // Residual
 Real refTemp = _refTempFunc.value(_t, _q_point[_qp]);

 return _transferCoeff * (_u[_qp] - refTemp) * _test[_i][_qp];
}

Real
NormalHeatFluxBC::computeQpJacobian()
{
 // Jacobian
 return _transferCoeff * _phi[_j][_qp] * _test[_i][_qp];
}
