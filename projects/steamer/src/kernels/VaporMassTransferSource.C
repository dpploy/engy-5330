//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VaporMassTransferSource.h"
/**
 * The forward declaration is so that we can declare the validParams() function
 * before we actually define the class... that way the definition isn't lost
 * at the bottom of the file.
 */

/// Forward Declarations
registerMooseObject("SteamerApp", VaporMassTransferSource);


/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template <> 
InputParameters validParams<VaporMassTransferSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($\ldots$), with the weak "
                             "form of $\ldots$.");
  params.addParam<Real>("sourceS",1.0,"Equation Term Coefficient");
  return params;
}

VaporMassTransferSource::VaporMassTransferSource(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _sourceS(getParam<Real>("sourceS"))
{
}

Real
VaporMassTransferSource::computeQpResidual()
{
  return  _sourceS *_test[_i][_qp] ;
}

Real
VaporMassTransferSource::computeQpJacobian()
{
  return 0.0;
}
