//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VaporConvection.h"

registerMooseObject("SteamerApp", VaporConvection);

template<>
InputParameters validParams<VaporConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Vapor convection term.");
  params.addParam<Real>("rhoV",1.0,"Vapor mass density");
  params.addParam<Real>("rhoL",1.0,"liquid mass density");
  params.addRequiredCoupledVar("vaporFraction", "Vapor fraction coupled variable.");
  return params;
}

VaporConvection::VaporConvection(const InputParameters & parameters): 
    Kernel(parameters),
	_rhoV(getParam<Real>("rhoV")),
	_rhoL(getParam<Real>("rhoL")),
	_vaporFraction(coupledValue("vaporFraction")),
	_gradVaporFraction(coupledGradient("vaporFraction"))
{
}

Real VaporConvection::computeQpResidual()
{
 Real alpha = _vaporFraction[_qp];
 Real alphaPrime = _gradVaporFraction[_qp](0);

 Real v = _u[_qp];
 Real vPrime = _grad_u[_qp](0);

 Real rho = alpha *_rhoV + (1-alpha)*_rhoL;

 Real theta = _test[_i][_qp];

 return _rhoV * (alphaPrime * v + alpha * vPrime) * theta;
}

Real
VaporConvection::computeQpJacobian()
{
 return 0.0;
}
