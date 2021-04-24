//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Convected.h"

registerMooseObject("FHRApp", Convected);

template<>
InputParameters validParams<Convected>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("density",1.0,"Fluid mass density");
  params.addRequiredParam<Real>("heatCapacity","Fluid heat capacity");
  params.addRequiredParam<Real>("velocity","Fluid velocity");
  return params;
}

Convected::Convected(const InputParameters & parameters):
    Kernel(parameters),
    density(getParam<Real>("density")),
    heatCapacity(getParam<Real>("heatCapacity")),
    velocity(getParam<Real>("velocity"))
{
}

Real
Convected::computeQpResidual()
{
  // Residual
  Real blin =  - density * heatCapacity * velocity * _grad_u[_qp] * _test[_i][_qp];
  //printf("%f ", _qp);
  return blin;
}

Real
Convected::computeQpJacobian()
{
  // Jacobian diagonal
  Real din = - density * heatCapacity * velocity * _grad_phi[_j][_qp] * _test[_i][_qp];
  //printf("%f ", _qp);
  return din;
}
