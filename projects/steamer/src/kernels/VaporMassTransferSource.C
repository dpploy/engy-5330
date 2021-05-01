//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VaporMassTransferSource.h"

registerMooseObject("SteamerApp", VaporMassTransferSource);

template <>
InputParameters validParams<VaporMassTransferSource>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Vapor mass transfer source.");
  params.addParam<Real>("sourceS",1.0,"Vapor mass transfer source");
  return params;
}

VaporMassTransferSource::VaporMassTransferSource(const InputParameters & parameters): 
    Kernel(parameters),
    _sourceS(getParam<Real>("sourceS"))
{
}

Real
VaporMassTransferSource::computeQpResidual()
{
  return - _sourceS *_test[_i][_qp] ;
}

Real
VaporMassTransferSource::computeQpJacobian()
{
  return 0.0;
}
