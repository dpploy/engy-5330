//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

/**
 * The forward declaration is so that we can declare the validParams() function
 * before we actually define the class... that way the definition isn't lost
 * at the bottom of the file.
 */

/// Forward Declarations
class VaporDriftDiffusion;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template <>
InputParameters validParams<VaporDriftDiffusion>();

/**
 * This kernel implements the following operator:
 *
 * $ u ... v $
 *
 * where v is a test function and u is an admissible solution
 */

/// VaporDriftDiffusion class inherits from Kernel class
class VaporDriftDiffusion : public Kernel
{
public:

  /**
   * This is the constructor declaration.  This class takes a
   * InputParameters object, just like other
   * Kernel-derived classes.
   */
  VaporDriftDiffusion(const InputParameters & parameters);

protected:
  /// Required residual for standard kernels in MOOSE
  virtual Real computeQpResidual() override;

  /// Required Jacobian for standard kernels in MOOSE
  /** This function returns the diagonal of the Jacobian to be used as a preconditioner
   * in the linear sub-problem.
   */
  virtual Real computeQpJacobian() override;

  /// The variables which holds the value for the VaporDriftDiffusion coefficient
  const Real _diffCoeff;
  const Real _rhoV;
  const Real _rhoL;
  const VariableValue & _fractionVapor;
  const VariableGradient & _grad_fractionVapor;
  //const VariableValue & _velocityMixture;
  //const VariableGradient & _grad_velocityMixture;
  /// The variables which holds the value for the MixtureMassBalDivergence coefficient

};
