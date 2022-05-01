//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SourceTerm.h"

registerMooseObject("NeutronBallApp", SourceTerm);

template<>
InputParameters validParams<SourceTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Source term kernel");
  params.addRequiredCoupledVar("coupledGroupA", "Coupled group A.");
  params.addRequiredCoupledVar("coupledGroupB", "Coupled group B.");
  params.addParam<Real>("sourceS",0.0,"Source term");
  params.addParam<Real>("sigma_sa",0.0,"Scatter to A");
  params.addParam<Real>("sigma_sb",0.0,"Scatter to B");
  return params;
}

SourceTerm::SourceTerm(const InputParameters & parameters): 
    Kernel(parameters),
    _coupledGroupA(coupledValue("coupledGroupA")),
    _coupledGroupB(coupledValue("coupledGroupB")),
    _sourceS(getParam<Real>("sourceS")),
    _sigma_sa(getParam<Real>("sigma_sa")),
    _sigma_sb(getParam<Real>("sigma_sb"))
{
}

Real
SourceTerm::computeQpResidual()
{
 return (_sigma_sa * _coupledGroupA[_qp] + _sigma_sb * _coupledGroupB[_qp]) * _test[_i][_qp];
}

Real
SourceTerm::computeQpJacobian()
{
 return (_sigma_sa + _sigma_sb) * _phi[_j][_qp] * _test[_i][_qp];
}
