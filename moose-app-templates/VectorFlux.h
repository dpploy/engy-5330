//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

class VectorFlux;

template <>
InputParameters validParams<VectorFlux>();

/**
 * Auxiliary kernel responsible for computing the flux of a variable.
 */
class VectorFlux : public VectorAuxKernel
{
public:
  static InputParameters validParams();

  VectorFlux(const InputParameters & parameters);

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual RealVectorValue computeValue() override;

  /// The gradient of a coupled variable
  const VariableGradient & _gradientVariableName;

  /// Add here other parameters needed
  Real _param1;
  Real _param2;
};  
