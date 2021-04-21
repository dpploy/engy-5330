//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MixtureMassBalDivergence.h"

/**
 * All MOOSE based object classes you create must be registered using this macro.
 * The first argument is the name of the App you entered in when running the stork.sh
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("SteamerApp", MixtureMassBalDivergence);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<MixtureMassBalDivergence>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($\ldots$), with the weak "
                             "form of $\ldots$.");
  params.addParam<Real>("rho_v",1.0,"Vapor density");
  params.addParam<Real>("rho_l",1.0,"liquid density");
  return params;
}

MixtureMassBalDivergence::MixtureMassBalDivergence(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _rho_v(getParam<Real>("rho_v")),
	_rho_l(getParam<Real>("rho_l")),
	_fractionVapor(coupledValue("fractionVapor")),
	_velocityMixture(coupledValue("velocityMixture")),
	_grad_fractionVapor(coupledGradient("fractionVapor")),
	_grad_velocityMixture(coupledGradient("velocityMixture"))
{
}

Real
MixtureMassBalDivergence::computeQpResidual()
{
  return (_grad_fractionVapor[_qp] *(_rho_v-_rho_l)*_velocityMixture[_qp]+(_fractionVapor[_qp]*_rho_v+(1-_fractionVapor[_qp])*_rho_l)*_grad_velocityMixture[_qp])*_test[_i][_qp];
}

Real
MixtureMassBalDivergence::computeQpJacobian()
{
	return 0.0;
}
