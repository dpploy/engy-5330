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

#pragma once

#include "SideIntegralPostprocessor.h"
#include "MooseVariableInterface.h"

class ScalarVariableSideIntegralPostprocessor : public SideIntegralPostprocessor,
                                          public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  ScalarVariableSideIntegralPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  /// Holds a parameter
  const Real _param1;
};
