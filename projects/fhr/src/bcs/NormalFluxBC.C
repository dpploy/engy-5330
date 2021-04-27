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
  InputParameters params = IntegratedBC::validParams();
  params.set<Real>("convectionCoeff") = 0;
  params.set<Real>("bias") = 0;
  params.set<Real>("zmax") = 0;
  params.addRequiredParam<FunctionName>("function", "the function");
  params.addClassDescription(
   "Parameters for the flux boundary condition. Rename and/or add as necessary.");
  return params;
}

NormalFluxBC::NormalFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _convectionCoeff(getParam<Real>("convectionCoeff")),
    _bias(getParam<Real>("bias")),
    _zmax(getParam<Real>("zmax")),
    _func(getFunction("function"))
{
}

Real
NormalFluxBC::computeQpResidual()
{
 // Implement the return:
 // e.g. return _param1 * (_u[_qp] - _param2) * _test[_i][_qp];
  Real funcValue = 600 + 100 * (_func.value(_t, _q_point[_qp]) / _zmax);
  return _convectionCoeff * (_u[_qp] - funcValue) * _test[_i][_qp];
  //FIXME return 0.0; // remove this line
}

Real
NormalFluxBC::computeQpJacobian()
{
 // Implement the return:
 // e.g. return _param1 * _phi[_j][_qp] * _test[_i][_qp];
 // FIXME return 0.0; // remove this line
 //Real funcValue = 600 + 100 * (_func.value(_t, _q_point[_qp]) / _zmax);
 return _convectionCoeff * _phi[_j][_qp] * _test[_i][_qp];
}
