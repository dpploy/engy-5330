//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "NormalFluxBC.h"

registerMooseObject("Engy5310P1App", NormalFluxBC);

defineLegacyParams(NormalFluxBC);

InputParameters
NormalFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addParam<Real>("transferCoeff","Mass Transfer Coefficient");
  params.addParam<Real>("reference","Ambient Concentration");
  params.addClassDescription(
   "Parameters for the flux boundary condition. Rename and/or add as necessary.");
  return params;
}

NormalFluxBC::NormalFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _transferCoeff(getParam<Real>("transferCoeff")),
    _reference(getParam<Real>("reference"))
{
}

Real
NormalFluxBC::computeQpResidual()
{
  return _transferCoeff * (_u[_qp] - _reference) * _test[_i][_qp];
  
}

Real
NormalFluxBC::computeQpJacobian()
{

 return _transferCoeff * _phi[_j][_qp] * _test[_i][_qp];
  
}