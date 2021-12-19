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

#pragma once

#include "ElementIntegralPostprocessor.h"

class ScalarVariableIntegralPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  ScalarVariableIntegralPostprocessor(const InputParameters & parameters);

  virtual Real getValue();

protected:
  virtual Real computeQpIntegral() override;

  Real _param1;
  Real _param2;

  const VariableValue & _variableName1;
  const VariableValue & _variableName2;
};
