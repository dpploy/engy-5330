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
class DiffusionFluxComponent;
/**
 * Auxiliary kernel responsible for computing the flux of a variable.
*/
template <>
InputParameters validParams<DiffusionFluxComponent>();

class DiffusionFluxComponent : public AuxKernel
{
public:
  static InputParameters validParams();

  DiffusionFluxComponent(const InputParameters & parameters);

protected:
  /**
   * AuxKernels MUST override computeValue.  computeValue() is called on
   * every quadrature point.  For Nodal Auxiliary variables those quadrature
   * points coincide with the nodes.
   */
  virtual Real computeValue() override;
  
  /// Will hold 0, 1, or 2 corresponding to x, y, or z.
  int _component;

  /// The derivative of a coupled variable
  const VariableGradient & _gradU;

  /// Add here other parameters needed
  Real _diffCoeff;
  //Real _grad_test;
};  