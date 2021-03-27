//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusionTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro.
 * The first argument is the name of the App you entered in when running the stork.sh
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("NeutronBallApp", DiffusionTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<DiffusionTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The scaled Laplacian operator (SD\\nabla \\cdot \\nabla u$) with the weak form of S(D\\nabla)");
  params.addParam<Real>("diffCoeff",0.0,"Equation Term Coefficient");
  return params;
}

DiffusionTerm::DiffusionTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _diffCoeff(getParam<Real>("diffCoeff"))
{
}

Real
DiffusionTerm::computeQpResidual()
{
  return - _diffCoeff * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
DiffusionTerm::computeQpJacobian()
{
  return - _diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
