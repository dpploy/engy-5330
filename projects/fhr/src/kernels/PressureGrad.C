//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PressureGrad.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("FHRApp", PressureGrad);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<PressureGrad>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("diffCoeff",1.0,"Equation Term Coefficient");
  return params;
}

PressureGrad::PressureGrad(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    PressureGrad(getParam<Real>("diffCoeff"))
{
}

Real
PressureGrad::computeQpResidual()
{
  // Implement the return
  Real blin =  -1 *  _u[_qp] * _grad_test[_i][_qp];
  //printf("%f ", _qp);
  return blin;
  //FIXME return 0.0 // remove this line
}

Real
PressureGrad::computeQpJacobian()
{
  // Implement the return
  Real din = -1 * _phi[_j][_qp] * _grad_test[_i][_qp];
  //printf("%f ", _qp);
  return din;

  //FIXME return 0.0 // remove this line
}
