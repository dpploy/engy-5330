//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SourceTerm.h"

registerMooseObject("FuelPlateApp", SourceTerm);

template<>
InputParameters validParams<SourceTerm>()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Weak form of Source term");
  params.addParam<Real>("sourceSlope",0.0,"Source slope term");
  return params;
}

SourceTerm::SourceTerm(const InputParameters & parameters) : 
    Kernel(parameters),
    _sourceSlope(getParam<Real>("sourceSlope"))
{
}

Real
SourceTerm::computeQpResidual()
{
 return - _sourceSlope * _u[_qp] * _test[_i][_qp];
}

Real
SourceTerm::computeQpJacobian()
{
 return - _sourceSlope * _phi[_j][_qp] * _test[_i][_qp];
}
