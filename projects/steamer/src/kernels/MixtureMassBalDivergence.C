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
  params.addParam<Real>("rhoV",1.0,"Vapor density");
  params.addParam<Real>("rhoL",1.0,"liquid density");
  params.addRequiredCoupledVar("velocityMixture", "Mixture velocity");
  return params;
}

MixtureMassBalDivergence::MixtureMassBalDivergence(const InputParameters & parameters):
    Kernel(parameters),
    _rhoV(getParam<Real>("rhoV")),
	_rhoL(getParam<Real>("rhoL")),
	_velocityMixture(coupledValue("velocityMixture")),
	_gradVelocityMixture(coupledGradient("velocityMixture"))

{
}

Real
MixtureMassBalDivergence::computeQpResidual()
{
 Real alpha = _u[_qp];
 Real alphaPrime = _grad_u[_qp](0);

 Real v = _velocityMixture[_qp];
 Real vPrime = _gradVelocityMixture[_qp](0);

 Real rho = alpha *_rhoV + (1-alpha)*_rhoL;
 Real delRho = _rhoV - _rhoL;

 Real theta = _test[_i][_qp];

 return (alphaPrime * delRho * v + rho * vPrime) * theta;
}

Real
MixtureMassBalDivergence::computeQpJacobian()
{
 //return ((_grad_phi[_j][_qp] *(_rhoV-_rhoL)*_velocityMixture[_qp])+ _gradVelocityMixture[_qp]*(_u[_qp]*_rhoV+_rhoL-_phi[_j][_qp]*_rhoL))*_test[_i][_qp];
 return 0.0;
}
