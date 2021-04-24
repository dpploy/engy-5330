//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "BoussinesqSource.h"

registerMooseObject("FHRApp", BoussinesqSource);

template<>
InputParameters validParams<BoussinesqSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("density",1.0,"fluid mass density");
  params.addRequiredParam<Real>("beta","fluid thermal expansion coefficient");
  params.addRequiredParam<Real>("refTemp","fluid thermal expansion reference temperature");
  return params;
}

BoussinesqSource::BoussinesqSource(const InputParameters & parameters): 
    Kernel(parameters),
    density(getParam<Real>("density")),
    beta(getParam<Real>("beta")),
    refTemp(getParam<Real>("refTemp")) 
{
}

// Residual
Real
BoussinesqSource::computeQpResidual()
{
  Real blin =  - (density + beta * refTemp) * _test[_i][_qp];
  //printf("%f ", _qp);
  return blin;
}

// Jacobian diagonal
Real
BoussinesqSource::computeQpJacobian()
{
  //Real din = - diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  //printf("%f ", _qp);
  //return din;

  return 0.0; // remove this line
}
