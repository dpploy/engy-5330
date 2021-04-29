//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ConvectionTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("Engy5310P1App", ConvectionTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<ConvectionTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<RealVectorValue>("velocity","Bulk Velocity");
  return params;
}

ConvectionTerm::ConvectionTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _velocity(getParam<RealVectorValue>("velocity"))
{
}

Real
ConvectionTerm::computeQpResidual()
{
 
  return  -_grad_u[_qp] *_velocity * _test[_i][_qp];

}

Real
ConvectionTerm::computeQpJacobian()
{

  return  -_grad_phi[_j][_qp] * _velocity * _test[_i][_qp];
  
}
