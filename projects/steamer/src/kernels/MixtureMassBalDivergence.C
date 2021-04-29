//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MixtureMassBalDivergence.h"

registerMooseObject("SteamerApp", MixtureMassBalDivergence);

template<>
InputParameters validParams<MixtureMassBalDivergence>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Divergence term of the mixture mass balance equation");
  params.addParam<Real>("rhoV",1.0,"Vapor mass density");
  params.addParam<Real>("rhoL",1.0,"liquid mass density");
  params.addRequiredCoupledVar("vaporFraction", "Vapor fraction");
  return params;
}

MixtureMassBalDivergence::MixtureMassBalDivergence(const InputParameters & parameters):
    Kernel(parameters),
    _rhoV(getParam<Real>("rhoV")),
    _rhoL(getParam<Real>("rhoL")),
    _vaporFraction(coupledValue("vaporFraction")),
    _gradVaporFraction(coupledGradient("vaporFraction"))
{
}

Real
MixtureMassBalDivergence::computeQpResidual()
{
 Real v = _u[_qp];
 Real vPrime = _grad_u[_qp](0);

 Real alpha = _vaporFraction[_qp];
 Real alphaPrime = _gradVaporFraction[_qp](0);

 Real rho = alpha * _rhoV + (1-alpha) * _rhoL;
 Real delRho = _rhoV - _rhoL;

 Real theta = _test[_i][_qp];

 return (delRho * alphaPrime * v + rho * vPrime) * theta;
}

Real
MixtureMassBalDivergence::computeQpJacobian()
{
 return 0.0;
}
