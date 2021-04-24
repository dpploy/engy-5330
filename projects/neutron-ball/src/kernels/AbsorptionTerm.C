//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AbsorptionTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro.
 * The first argument is the name of the App you entered in when running the stork.sh
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("NeutronBallApp", AbsorptionTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<AbsorptionTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The weak form of the neutron intensity absoption term");
  params.addParam<Real>("sigmaA",0.0,"Macroscopic absoption cross-section");
  return params;
}

AbsorptionTerm::AbsorptionTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _sigmaA(getParam<Real>("sigmaA"))
{
}

Real
AbsorptionTerm::computeQpResidual()
{
  return - ( - _sigmaA * _u[_qp] * _test[_i][_qp] );
  //return _sigmaA * _test[_i][_qp];
}

Real
AbsorptionTerm::computeQpJacobian()
{
  return - ( - _sigmaA * _phi[_j][_qp] * _test[_i][_qp] );
  //return 0.0;
}
