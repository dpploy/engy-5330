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

#include "SourceTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("FIRESBrickApp", SourceTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<SourceTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("sourceS",1.0,"Equation Term Coefficient");
  return params;
}

SourceTerm::SourceTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _sourceS(getParam<Real>("sourceS"))
{
}

Real
SourceTerm::computeQpResidual()
{
  
return _sourceS * _test[_i][_qp];

}

Real
SourceTerm::computeQpJacobian()
{
  return 0.0;
}
