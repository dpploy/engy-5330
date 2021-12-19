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
//
// Engy-5310: Computational Continuum Transport Phenomena
// UMass Lowell, Nuclear Chemical Engineering
// https://github.com/dpploy/engy-5310

#include "FluxComponent.h"

registerMooseObject("FIRESBrickApp", FluxComponent);

defineLegacyParams(FluxComponent);

InputParameters
FluxComponent::validParams()
{
  InputParameters params = AuxKernel::validParams();
  MooseEnum component("x y z");
  params.addClassDescription("Compute one component of the flux of a field.");
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of flux.");
  // Add a "coupling parameter" to get a variable from the input file.
  params.addRequiredCoupledVar("field", "The variable field.");
  // Add add other parameters from input file
  params.addRequiredParam<Real>("param1", "Parameter 1 meaning");
  params.addRequiredParam<Real>("param2", "Parameter 2 meaning");

  return params;
}

FluxComponent::FluxComponent(const InputParameters & parameters):
    AuxKernel(parameters),
    _component(getParam<MooseEnum>("component")),
    // Initialize variable gradient
    _gradVariableName_component(coupledGradient("field")),
    // Initialize parameters
    _param1(getParam<Real>("param1")),
    _param2(getParam<Real>("param2"))
{
}

Real
FluxComponent::computeValue()
{
  // Access the gradient of the variable at this quadrature point, then pull out the 
  // "component" of it requested (x, y or z). Note, that getting a particular component 
  // of a gradient is done using the parenthesis operator.
  return - _param1 * _param2 * _gradVariableName_component[_qp](_component);
}
