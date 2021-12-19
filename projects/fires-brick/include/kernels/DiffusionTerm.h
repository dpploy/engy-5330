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

#include "Kernel.h"

class DiffusionTerm : public Kernel
{
public:

  DiffusionTerm(const InputParameters & parameters);

protected:
  /// Residual
  virtual Real computeQpResidual() override;

  /// Jacobian diagonal
  virtual Real computeQpJacobian() override;

  /// User variables e.g.
  const Real _diffCoeff;
//const Real alfa;
//const Real beta;

};
