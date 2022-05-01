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

#include "HeatConduction.h"

registerMooseObject("HeatedPlateApp", HeatConduction);

template<>
InputParameters validParams<HeatConduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Weak form of heat conduction");
  params.addParam<Real>("thermCond", 1.0, "Thermal conductivity coefficient");
  return params;
}

HeatConduction::HeatConduction(const InputParameters & parameters): 
  Kernel(parameters),
  _thermCond(getParam<Real>("thermCond"))
{
}

Real
HeatConduction::computeQpResidual()
{
  RealVectorValue q = - _thermCond * _grad_u[_qp];

 // Residual
  return - q * _grad_test[_i][_qp];
}

Real
HeatConduction::computeQpJacobian()
{
  // Jacobian diagonal
  return _thermCond * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
