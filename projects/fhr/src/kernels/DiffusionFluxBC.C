//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusionFluxBC.h"

registerMooseObject("FHRApp", DiffusionFluxBC);

defineLegacyParams(DiffusionFluxBC);

InputParameters
DiffusionFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.set<Real>("conductivity") = 0;
  params.set<Real>("convectivity") = 0;
  params.set<Real>("refTemp") = 0;
  params.addClassDescription(
   "Parameters for the flux boundary condition. Rename and/or add as necessary.");
  return params;
}

DiffusionFluxBC::DiffusionFluxBC(const InputParameters & parameters): 
    IntegratedBC(parameters),
    _conductivity(getParam<Real>("conductivity")),
    _convectivity(getParam<Real>("convectivity")),
    _refTemp(getParam<Real>("refTemp"))
{
}

Real
DiffusionFluxBC::computeQpResidual()
{
 // Implement the return:
 return (conductivity / convectivity) * (_u[_qp] - refTemp) * _test[_i][_qp];
 // FIXME return 0.0; // remove this line
}

Real
EquationFluxBC::computeQpJacobian()
{
 // Implement the return:
 return (conductivity / convectivity) * _phi[_j][_qp] * _test[_i][_qp];
  //FIXME return 0.0; // remove this line
}
