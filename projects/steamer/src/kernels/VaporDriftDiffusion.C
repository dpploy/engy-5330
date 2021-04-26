//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VaporDriftDiffusion.h"

/**
 * All MOOSE based object classes you create must be registered using this macro.
 * The first argument is the name of the App you entered in when running the stork.sh
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("SteamerApp", VaporDriftDiffusion);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<VaporDriftDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak "
                             "form of $...$.");
  params.addParam<Real>("diffCoeff",1.0,"Equation Term Coefficient");
  params.addParam<Real>("rhoV",1.0,"Vapor density");
  params.addParam<Real>("rhoL",1.0,"liquid density");
  params.addRequiredCoupledVar("some_variable", "The gradient of this variable will be used as the velocity vector.");
  return params;
}

VaporDriftDiffusion::VaporDriftDiffusion(const InputParameters & parameters): 
    Kernel(parameters),
    // Set the coefficient for the equation term

  _diffCoeff(getParam<Real>("diffCoeff")),
	_rhoV(getParam<Real>("rhoV")),
	_rhoL(getParam<Real>("rhoL")),
	_fractionVapor(coupledValue("some_variable")),
	//_velocityMixture(coupledValue("velocityMixture")),
	_grad_fractionVapor(coupledGradient("some_variable"))
	//_grad_velocityMixture(coupledGradient("velocityMixture"))

{
}

Real
VaporDriftDiffusion::computeQpResidual()
{
	return ((_grad_fractionVapor[_qp](0)*_rhoV*_u[_qp] + _fractionVapor[_qp]*_rhoV*_grad_u[_qp](0)) - ((_rhoL*_fractionVapor[_qp]*_diffCoeff*_rhoV*_grad_fractionVapor[_qp](0)*(_rhoV-_rhoL))/((_fractionVapor[_qp]*_rhoV+(1-_fractionVapor[_qp])*_rhoL)*(_fractionVapor[_qp]*_rhoV+(1-_fractionVapor[_qp])*_rhoL)))) *_test[_i][_qp];
}

Real
VaporDriftDiffusion::computeQpJacobian()
{
	return 0.0;
}
