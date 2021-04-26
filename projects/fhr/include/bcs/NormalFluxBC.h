//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADFunctionNeumannBC.h"

// Forward Declarations
class NormalFluxBC;

template <>
InputParameters validParams<NormalFluxBC>();

class NormalFluxBC : public ADFunctionNeumannBC
{
public:
  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  static InputParameters validParams();

  NormalFluxBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  virtual ADReal computeQpJacobian();

private:
  /**
   * Parameters of the flux boundary condition. Rename and/or add parameter names below.
   */
  Real _conductivity;
  Real _convectivity;
  Real _zmax;
  const Function & _func;
};
