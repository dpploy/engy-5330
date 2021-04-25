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

class NuclearHeat : public Kernel
{
public:

  NuclearHeat(const InputParameters & parameters);

  // Residual
  virtual Real computeQpResidual() override;

  // Jacobian diagonal
  virtual Real computeQpJacobian() override;

  // Nuclear heating power
  Real sourceS;
};
