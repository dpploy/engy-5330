//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VaporDriftFluxDiffusion.h"

registerMooseObject("SteamerApp", VaporDriftFluxDiffusion);

template<>
InputParameters validParams<VaporDriftFluxDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Diffusion drift flux term.");
  params.addParam<Real>("diffCoeff",1.0,"Drift flux diffusion coefficient");
  params.addParam<Real>("rhoV",1.0,"Vapor mass density");
  params.addParam<Real>("rhoL",1.0,"liquid mass density");
  params.addRequiredCoupledVar("velocity", "Mixture velocity coupled variable.");
  return params;
}

VaporDriftFluxDiffusion::VaporDriftFluxDiffusion(const InputParameters & parameters): 
    Kernel(parameters),
    _diffCoeff(getParam<Real>("diffCoeff")),
	_rhoV(getParam<Real>("rhoV")),
	_rhoL(getParam<Real>("rhoL")),
	_velocity(coupledValue("velocity")),
	_gradVelocity(coupledGradient("velocity"))
{
}

Real VaporDriftFluxDiffusion::computeQpResidual()
{
 Real alpha = _u[_qp];
 Real alphaPrime = _grad_u[_qp](0);

 Real v = _velocity[_qp];
 Real vPrime = _gradVelocity[_qp](0);

 Real rho = alpha *_rhoV + (1-alpha)*_rhoL;
 Real delRho = _rhoV - _rhoL;

 Real wL = _rhoL/rho;
 Real wV = _rhoV/rho;

 Real theta = _test[_i][_qp];

 return - alpha * wL * wV * delRho * _diffCoeff * alphaPrime * theta;
}

Real
VaporDriftFluxDiffusion::computeQpJacobian()
{
 return 0.0;
}
