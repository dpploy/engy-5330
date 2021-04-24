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
 * This kernel implements the scaled Laplacian operator:
 *
 * $\nabla u \cdot \nabla v$
 *
 * where v is test function and u is an admissible solution
 */

/// DiffusionTerm class inherits from Kernel class
class DiffusionTerm : public Kernel
{
public:
  DiffusionTerm(const InputParameters & parameters);

protected:
  // Diffusion residual
  virtual Real computeQpResidual() override;

  // Diffusion Jacobian diagonal
  virtual Real computeQpJacobian() override;

  // Diffusion data member
  const Real _diffCoeff;
};
