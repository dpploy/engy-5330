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

#include "DiffusionTerm.h"

registerMooseObject("FuelPlateApp", DiffusionTerm);

template<>
InputParameters validParams<DiffusionTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Weak form of neutron diffusion");
  params.addParam<Real>("diffCoeff", 1.0, "Diffusion coefficient");
  return params;
}

DiffusionTerm::DiffusionTerm(const InputParameters & parameters): 
  Kernel(parameters),
  _diffCoeff(getParam<Real>("diffCoeff"))
{
}

Real
DiffusionTerm::computeQpResidual()
{
  RealVectorValue q = - _diffCoeff * _grad_u[_qp];

 // Residual
  return - q * _grad_test[_i][_qp];
}

Real
DiffusionTerm::computeQpJacobian()
{
  // Jacobian diagonal
  return _diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
