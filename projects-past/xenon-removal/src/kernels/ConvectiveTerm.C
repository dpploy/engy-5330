//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConvectiveTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("Engy5310p2App", ConvectiveTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<ConvectiveTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<RealVectorValue>("velocity", "Velocity Vector");
  params.addParam<Real>("voidFraction",1.0,"DiffusionTerm Coefficient");

  return params;
}

ConvectiveTerm::ConvectiveTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
     _velocity(getParam<RealVectorValue>("velocity")),
	 _voidFraction(getParam<Real>("voidFraction"))

{
}

Real
ConvectiveTerm::computeQpResidual()
{
 return  _voidFraction*_grad_u[_qp] * _velocity * _test[_i][_qp];
}

Real
ConvectiveTerm::computeQpJacobian()
{
return  _voidFraction* _grad_phi[_j][_qp] * _velocity * _test[_i][_qp];
}