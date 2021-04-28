//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DiffusionFluxComponent.h"

registerMooseObject("Engy5310P1App", DiffusionFluxComponent);

defineLegacyParams(DiffusionFluxComponent);

InputParameters
DiffusionFluxComponent::validParams()
{
  InputParameters params = AuxKernel::validParams();
  MooseEnum component("x y z");
  params.addClassDescription("Compute one component of the flux of a field.");
  params.addRequiredParam<MooseEnum>("component", component, "The desired component of flux.");
  // Add a "coupling parameter" to get a variable from the input file.
  params.addRequiredCoupledVar("field", "The variable field.");
  // Add add other parameters from input file
  params.addRequiredParam<Real>("diffCoeff", "Parameter 1 meaning");
  //params.addRequiredParam<Real>("grad_test", "Parameter 2 meaning");

  return params;
}

DiffusionFluxComponent::DiffusionFluxComponent(const InputParameters & parameters)
  : AuxKernel(parameters),
    _component(getParam<MooseEnum>("component")),
    // Initialize variable gradient
    _gradU(coupledGradient("field")),
    // Initialize parameters
    _diffCoeff(getParam<Real>("diffCoeff"))
    //_grad_test(getParam<Real>("grad_test"))
{
}

Real
DiffusionFluxComponent::computeValue()
{
  // Access the gradient of the variable at this quadrature point, then pull out the 
  // "component" of it requested (x, y or z). Note, that getting a particular component 
  // of a gradient is done using the parenthesis operator.
 return - _diffCoeff * _gradU[_qp](_component);
}