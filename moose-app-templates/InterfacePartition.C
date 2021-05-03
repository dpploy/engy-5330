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
//
// https://github.com/dpploy/engy-5310

#include "InterfacePartition.h"

registerMooseObject("Engy5310App-FIXME", InterfacePartition);

defineLegacyParams(InterfacePartition);

InputParameters
InterfacePartition::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("kCoeff", "Partition coefficient k >= 0");
  params.addClassDescription("Partition u_primary = k u_neighbor");
  return params;
}

InterfacePartition::InterfacePartition(const InputParameters & parameters): 
    InterfaceKernel(parameters), 
    _kCoeff(getParam<Real>("kCoeff"))
{
}

Real
InterfacePartition::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    // Primary residual = u_primary - k * u_neighbor
    // Weak form for primary domain is: (u_primary - k*u_neighbor, test)
    case Moose::Element:
      r = _test[_i][_qp] * (_u[_qp] - _kCoeff * _neighbor_value[_qp]);
      break;

    // Secondary residual: -(u_primary - k*u_neighbor, test),
    // flip the sign because the integration direction is opposite.
    case Moose::Neighbor:
      r = - _test_neighbor[_i][_qp] * (_u[_qp] - _kCoeff * _neighbor_value[_qp]);
      break;
  }
  return r;
}

Real
InterfacePartition::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::NeighborNeighbor:
      jac = -_test_neighbor[_i][_qp] * -_kCoeff * _phi_neighbor[_j][_qp];
      break;
    case Moose::NeighborElement:
      jac = -_test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * -_kCoeff * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}