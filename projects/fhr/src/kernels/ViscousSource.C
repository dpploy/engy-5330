//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ViscousSource.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("Engy5310p1App", ViscousSource);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<ViscousSouce>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("viscocity",1.0,"Equation Term Coefficient");
  params.addParam<Real>("A");
  params.addPaam<Real>("B");
  return params;
}

ViscousSource::ViscousSource(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    viscocity(getParam<Real>("viscocity")),
    A(getParam<Real>("A")),
    B(getParam<Real>("B"));
{
}

Real
ViscousSource::computeQpResidual()
{
  // Implement the return
  Real blin =  - viscocity * (4 * A^2 * x + 4 * A * B) * _test[_i][_qp];
  //printf("%f ", _qp);
  return blin;
  //FIXME return 0.0 // remove this line
}

Real
ViscousSource::computeQpJacobian()
{
  // Implement the return
  //Real din = - diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  //printf("%f ", _qp);
  //return din;

  return 0.0; // remove this line
}
