

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
class SourceTerm;

/**
 * validParams returns the parameters that this Kernel accepts / needs
 * The actual body of the function MUST be in the .C file.
 */
template <>
InputParameters validParams<SourceTerm>();

/**
 * This kernel implements the Laplacian operator:
 *
 * $\nabla u \cdot \nabla v$
 *
 * where v is test function and u is an admissible solution
 */

/// EquationTerm class inherits from Kernel class
class SourceTerm : public Kernel
{
public:

  /**
   * This is the constructor declaration.  This class takes a
   * InputParameters object, just like other
   * Kernel-derived classes.
   */
  SourceTerm(const InputParameters & parameters);

protected:
  /// Required residual for standard kernels in MOOSE
  virtual Real computeQpResidual() override;

  /// Required Jacobian for standard kernels in MOOSE
  /** This function returns the diagonal of the Jacobian to be used as a preconditioner
   * in the linear sub-problem.
   */
  virtual Real computeQpJacobian() override;

  /// The variables which holds the value for the EquationTerm coefficient
  const Real _massArea;
  const Real _henryC;
  const Real _massC;
  const Real _heC;
  const VariableValue & _uCoupled;
  
};