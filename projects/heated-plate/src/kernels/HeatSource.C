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

#include "HeatSource.h"

registerMooseObject("HeatedPlateApp", HeatSource);

template<>
InputParameters validParams<HeatSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Heat source term");
  params.addParam<Real>("sourceS", 1.0, "Heating power density");
  return params;
}

HeatSource::HeatSource(const InputParameters & parameters): 
  Kernel(parameters),
  _sourceS(getParam<Real>("sourceS"))
{
}

Real
HeatSource::computeQpResidual()
{
  // Residual
  return - _sourceS * _test[_i][_qp];
}

Real
HeatSource::computeQpJacobian()
{
  // Jacobian diagonal
  return 0.0;
}
