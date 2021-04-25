//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatConduction.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 */
registerMooseObject("FHRApp", HeatConduction);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<HeatConduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Weak form of thermal conductivity");
  params.addParam<Real>("thermCond",1.0,"Thermal conductivity coefficient");
  return params;
}

HeatConduction::HeatConduction(const InputParameters & parameters): 
    Kernel(parameters),
    _thermCond(getParam<Real>("thermCond"))
{
}

Real
HeatConduction::computeQpResidual()
{
  return _thermCond * _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
HeatConduction::computeQpJacobian()
{
  return _thermCond * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
