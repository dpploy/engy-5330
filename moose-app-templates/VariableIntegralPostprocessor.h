//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"

// Forward Declarations
// none

class VariableIntegralPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  VariableIntegralPostprocessor(const InputParameters & parameters);
  virtual Real getValue();

protected:
  virtual Real computeQpIntegral();

  Real _param1;
  Real _param2;

  const VariableValue & _variableName1;
  const VariableValue & _variableName2;
};
