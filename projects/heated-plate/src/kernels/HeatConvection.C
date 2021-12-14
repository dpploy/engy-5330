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

#include "HeatConvection.h"

registerMooseObject("HeatedPlateApp", HeatConvection);

template<>
InputParameters validParams<HeatConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Heat convection term.");
  params.addParam<Real>("massDensity", 1.0, "Fluid mass density");
  params.addParam<Real>("heatCapacity", 1.0, "Fluid heat capacity");
  params.addRequiredParam<RealVectorValue>("velocity", "Fluid velocity");
  return params;
}

HeatConvection::HeatConvection(const InputParameters & parameters):
  Kernel(parameters),
  _massDensity(getParam<Real>("massDensity")),
  _heatCapacity(getParam<Real>("heatCapacity")),
  _velocity(getParam<RealVectorValue>("velocity"))
{
}

Real
HeatConvection::computeQpResidual()
{
  // Residual
  return _massDensity * _heatCapacity * _grad_u[_qp] * _velocity * _test[_i][_qp];
}

Real
HeatConvection::computeQpJacobian()
{
  // Jacobian diagonal
  return _massDensity * _heatCapacity * _grad_phi[_j][_qp] * _velocity * _test[_i][_qp];
}
