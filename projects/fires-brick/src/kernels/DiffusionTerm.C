//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//
// Engy-5310: Computational Continuum Transport Phenomena
// UMass Lowell, Nuclear Chemical Engineering
// https://github.com/dpploy/engy-5310

#include "DiffusionTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("FIRESBrickApp", DiffusionTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<DiffusionTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("diffCoeff",1.0,"Equation Term Coefficient");
  params.addParam<Real>("alfa",1.0,"alfa value");
  params.addParam<Real>("beta",1.0,"beta value");
//  params.addRequiredParam<FunctionName>("refxFunc", "Profile function");
  
return params;
}

DiffusionTerm::DiffusionTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _diffCoeff(getParam<Real>("diffCoeff"))
  //alfa(getParam<Real>("alfa")),
// _refxFunc(getFunction("refxFunc")),
  //beta(getParam<Real>("beta"))

{
}

Real
DiffusionTerm::computeQpResidual()
{ 
//Real x = _refxFunc.value(_t, _q_point[_qp]);  

// Implement the return
    return - _diffCoeff * _grad_u[_qp] * _grad_test[_i][_qp];
//return - _alfa*_u[_qp] * _grad_u[_qp] * _grad_test[_i][_qp];
 // return -( beta+ alfa  * x ) * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
DiffusionTerm::computeQpJacobian()
{
  // Implement the return
  return - _diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
//return - _alfa*_phi[_j][_qp] * (_grad_phi[_j][_qp]+_alfa*_u*_grad_phi[_j][_qp]) * _grad_test[_i][_qp];
//return -( beta+ alfa * x ) * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
 }

