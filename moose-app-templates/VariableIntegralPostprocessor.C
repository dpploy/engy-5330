//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// Template includes
#include "VariableIntegralPostprocessor.h"

registerMooseObject("Engy5310PXApp-FIXME", VariableIntegralPostprocessor);

InputParameters
VariableIntegralPostprocessor::validParams()
{
  InputParameters params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("Computes entropy error.");
  params.addRequiredParam<Real>("param1", "Parameter 1 meaning");
  params.addRequiredParam<Real>("param2", "Parameter 2 meaning");
  params.addRequiredCoupledVar("variableName1", "meaning of variable 1");
  params.addRequiredCoupledVar("variableName2", "meaning of variable 2");

  return params;
}

VariableIntegralPostprocessor::VariableIntegralPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _param1(getParam<Real>("param1")),
    _param2(getParam<Real>("param2")),
    _variableName1(coupledValue("variableName1")),
    _variableName2(coupledValue("variableName2"))
{
}

Real
VariableIntegralPostprocessor::getValue()
{
  return ElementIntegralPostprocessor::getValue();
}

Real
VariableIntegralPostprocessor::computeQpIntegral()
{
 // e.g. Real integrand = (_variableName1[_qp] / _param1) * std::pow(_param1 / _variableName2[_qp], 1.4) - 1.;
 FIXME return integrand;
}
