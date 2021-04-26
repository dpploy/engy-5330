//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NormalFluxBC.h"
#include "Function.h"

registerMooseObject("Engy5310p1App", NormalFluxBC);

defineLegacyParams(NormalFluxBC);

InputParameters
NormalFluxBC::validParams()
{
  InputParameters params = ADFunctionNeumannBC::validParams();
  params.addRequiredParam<Real>("conductivity") = 0;
  params.addRequiredParam<Real>("convectivity") = 0;
  params.addRequiredParam<FunctionName>("function") = 0;
  params.addRequiredParam<Real>("zmax") = 0;
  params.addClassDescription(
   "Parameters for the flux boundary condition. Rename and/or add as necessary.");
  return params;
}

NormalFluxBC::NormalFluxBC(const InputParameters & parameters)
  : ADFunctionNeumannBC(parameters),
    _conductivity(getParam<Real>("conductivity")),
    _convectivity(getParam<Real>("convectivity")),
    _func(getFunction("function")),
    _zmax(getParam<Real>("zmax"))
{
}

ADReal
NormalFluxBC::computeQpResidual()
{
 // Implement the return:
 // help
 Real val = _func.value(_t, _q_point[_qp]);
 Real refTemp = 600 + (100 * (val/_zmax));
 return (- _conductivity / _convectivity) * (_u[_qp] - refTemp) * _test[_i][_qp];
}

ADReal
NormalFluxBC::computeQpJacobian()
{
 // Implement the return:
 return (- _conductivity / _convectivity) * _phi[_j][_qp] * _test[_i][_qp];
}
