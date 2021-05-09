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

#include "AuxKernel.h"

class FluxComponent : public AuxKernel
{
public:
  static InputParameters validParams();

  FluxComponent(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  
  /// Will hold 0, 1, or 2 component corresponding to x, y, or z.
  int _component;

  /// The derivative of a coupled variable
  const VariableGradient & _gradVariableName_component;

  /// Add here other parameters needed
  Real _param1;
  Real _param2;
};  
