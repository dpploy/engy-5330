//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Engy-5310: Computational Continuum Transport Phenomena
//* UMass Lowell, Nuclear Chemical Engineering
//* https://github.com/dpploy/engy-5310

#pragma once

#include "IntegratedBC.h"

class Function;

class NormalFluxBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  NormalFluxBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

private:

  // User variables
  Real _transferCoeff;
};
