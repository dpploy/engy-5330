//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Convected.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("Engy5310p1App", Convected);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<Convected>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("density",1.0,"Equation Term Coefficient");
  params.addParam<Real>("heatCapacity");
  params.addParam<Real>("velocity");
  return params;
}

Convected::Convected(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    Convected(getParam<Real>("density"));
    Convected(getParam<Real>("heatCapacity"));
    Convected(getParam<Real>("velocity"));

{
}

Real
Convected::computeQpResidual()
{
  // Implement the return
  Real blin =  - density * heatCapacity * velocity * _grad_u[_qp] * _test[_i][_qp];
  //printf("%f ", _qp);
  return blin;
  //FIXME return 0.0 // remove this line
}

Real
Convected::computeQpJacobian()
{
  // Implement the return
  Real din = - density * heatCapacity * velocity * _grad_phi[_j][_qp] * _test[_i][_qp];
  //printf("%f ", _qp);
  return din;

  //FIXME return 0.0 // remove this line
}
