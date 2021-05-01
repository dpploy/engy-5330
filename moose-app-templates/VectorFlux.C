//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//
// https://github.com/dpploy/engy-5310

#include "VectorFlux.h"

registerMooseObject("Engy5310App-FIXME", VectorFlux);

InputParameters
VectorFlux::validParams()
{
  InputParameters params = VectorAuxKernel::validParams();
  params.addClassDescription("Computes the flux of a variable.");
  // Add a "coupling parameter" to get a variable from the input file.
  params.addRequiredCoupledVar("variableName", "The variable field.");
  // Add add other parameters from input file
  params.addRequiredParam<Real>("param1", "Parameter 1 meaning");
  params.addRequiredParam<Real>("param2", "Parameter 2 meaning");

  return params;
}

VectorFlux::VectorFlux(const InputParameters & parameters)
  : VectorAuxKernel(parameters),
    // Initialize variable gradient
    _gradVariableName(coupledGradient("variableName")),
    // Initialize parameters
    _param1(getParam<Real>("param1")),
    _param2(getParam<Real>("param2"))
{
}

RealVectorValue
VectorFlux::computeValue()
{
  // Access the gradient of the variable at this quadrature point, then pull out the 
  // "component" of it requested (x, y or z). Note, that getting a particular component 
  // of a gradient is done using the parenthesis operator.
  FIXME return - _param1 * _param2 * _gradVariableName[_qp];
}
