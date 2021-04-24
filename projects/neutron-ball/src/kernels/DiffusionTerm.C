//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusionTerm.h"

registerMooseObject("NeutronBallApp", DiffusionTerm);

template<>
InputParameters validParams<DiffusionTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Weak form of the Neutron diffusion flux magnitude");
  params.addParam<Real>("diffCoeff",0.0,"Neutron diffusion coefficient");
  return params;
}

DiffusionTerm::DiffusionTerm(const InputParameters & parameters): 
    Kernel(parameters),
    _diffCoeff(getParam<Real>("diffCoeff"))
{
}

// Residual
Real
DiffusionTerm::computeQpResidual()
{
  return - (- _diffCoeff * _grad_u[_qp] * _grad_test[_i][_qp]);
}

// Jacobian diagonal
Real
DiffusionTerm::computeQpJacobian()
{
  return - (- _diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp]) ;
}
