//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//
// Engy-5310: Computational Continuum Transport Phenomena
// UMass Lowell, Nuclear Chemical Engineering
// https://github.com/dpploy/engy-5310

#include "EquationFluxBC.h"

registerMooseObject("FIRESBrickApp", EquationFluxBC);

defineLegacyParams(EquationFluxBC);

InputParameters
EquationFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.set<Real>("param1") = 0;
  params.set<Real>("param2") = 0;
  params.addClassDescription(
   "Parameters for the flux boundary condition. Rename and/or add as necessary.");
  return params;
}

EquationFluxBC::EquationFluxBC(const InputParameters & parameters):
    IntegratedBC(parameters),
    _param1(getParam<Real>("param1")),
    _param2(getParam<Real>("param2"))
{
}

Real
EquationFluxBC::computeQpResidual()
{
 // Implement the return:
  return _param1 * (_u[_qp] - _param2) * _test[_i][_qp];
}

Real
EquationFluxBC::computeQpJacobian()
{
 // Implement the return:
   return _param1 * _phi[_j][_qp] * _test[_i][_qp];
}
