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

class VectorFlux : public VectorAuxKernel
{
public:
  static InputParameters validParams();

  VectorFlux(const InputParameters & parameters);

protected:
  virtual RealVectorValue computeValue() override;

  /// The gradient of a coupled variable
  const VariableGradient & _gradientVariableName;

  /// Add here other parameters needed
  Real _param1;
  Real _param2;
};  
