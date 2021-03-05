//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ScalarVariableSideIntegralPostProcessor.h"

registerMooseObject("Engy5310PXApp-FIXME", ScalarVariableSideIntegralPostProcessor);

defineLegacyParams(ScalarVariableSideIntegralPostProcessor);

InputParameters
ScalarVariableSideIntegralPostProcessor::validParams()
{
  InputParameters params = SideIntegralPostprocessor::validParams();
  params.addClassDescription("Computes a surface integral of the specified variable");
  params.addRequiredCoupledVar("variable", "The name of the variable that this boundary condition applies to");
  params.addRequiredParam<Real>("param1", "Parameter 1 meaning");

  return params;
}

ScalarVariableSideIntegralPostProcessor::ScalarVariableSideIntegralPostProcessor(
    const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_ANY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _u(coupledValue("variable")),
    _grad_u(coupledGradient("variable")),
    _param1(getParam<Real>("param1"))
{
  addMooseVariableDependency(&mooseVariableField());
}

Real
ScalarVariableSideIntegralPostProcessor::computeQpIntegral()
{
 // e.g. return - _param1 * _grad_u[_qp] * _normals[_qp];
 FIXME return _u[_qp];
}
