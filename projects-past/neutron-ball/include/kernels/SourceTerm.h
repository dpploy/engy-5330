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

class SourceTerm : public Kernel
{
public:
  SourceTerm(const InputParameters & parameters);

protected:
  /// Residual
  virtual Real computeQpResidual() override;

  /// Jacobian diagonal
  virtual Real computeQpJacobian() override;

  /// User variables
  const VariableValue & _coupledGroupA;
  const VariableValue & _coupledGroupB;
  const Real _sourceS;
  const Real _sigma_sa;
  const Real _sigma_sb;

};
