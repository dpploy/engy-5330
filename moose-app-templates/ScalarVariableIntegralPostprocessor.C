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

// Template includes
#include "ScalarVariableIntegralPostprocessor.h"

registerMooseObject("Engy5310PXApp-FIXME", ScalarVariableIntegralPostprocessor);

InputParameters
ScalarVariableIntegralPostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Computes a volume integral of a scalar variable.");
  params.addRequiredParam<Real>("param1", "Parameter 1 meaning");
  params.addRequiredParam<Real>("param2", "Parameter 2 meaning");
  // Add a "coupling variable" to get a variable from the input file.
  params.addRequiredCoupledVar("variable1", "Variable to be used");
  params.addRequiredCoupledVar("variable2", "Variable to be used");

  return params;
}

ScalarVariableIntegralPostprocessor::ScalarVariableIntegralPostprocessor(const InputParameters & parameters):
    ElementIntegralPostprocessor(parameters),
    _param1(getParam<Real>("param1")),
    _param2(getParam<Real>("param2")),
    _variableName1(coupledValue("variable1")),
    _variableName2(coupledValue("variable2"))
{
}

Real
ScalarVariableIntegralPostprocessor::getValue()
{
  return ElementIntegralPostprocessor::getValue();
}

Real
ScalarVariableIntegralPostprocessor::computeQpIntegral()
{
 // e.g. Real integrand = (_variableName1[_qp] / _param1) * std::pow(_param1 / _variableName2[_qp], 1.4) - 1.;
 FIXME return integrand;
}
