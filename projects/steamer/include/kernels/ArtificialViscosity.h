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

class ArtificialViscosity : public Kernel
{
public:

  ArtificialViscosity(const InputParameters & parameters);

protected:
  // Residual
  virtual Real computeQpResidual() override;

  // Jacobian diagonal
  virtual Real computeQpJacobian() override;

  // Property variables
  const Real _viscosity;
};
