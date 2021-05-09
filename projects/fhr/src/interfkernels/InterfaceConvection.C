//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//
// Engy-5310: Computational Continuum Transport Phenomena
// UMass Lowell, Nuclear Chemical Engineering
// https://github.com/dpploy/engy-5310

#include "InterfaceConvection.h"

registerMooseObject("FHRApp", InterfaceConvection);

defineLegacyParams(InterfaceConvection);

InputParameters
InterfaceConvection::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addRequiredParam<Real>("convectionCoeff", "Partition coefficient k >= 0");
  params.addRequiredParam<Real>("conductivity", "primary conductivity");
  params.addRequiredParam<Real>("neighborConductivity", "conductivity in neighboring phase");
  params.addClassDescription("h(u-u_neighbor) = k  /del  u");
  return params;
}

InterfaceConvection::InterfaceConvection(const InputParameters & parameters): 
    InterfaceKernel(parameters), 
    _convectionCoeff(getParam<Real>("convectionCoeff")),
    _conductivity(getParam<Real>("conductivity")),
    _neighborConductivity(getParam<Real>("neighborConductivity"))
{
}

Real
InterfaceConvection::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    // Primary residual = h * (u_primary - u_neighbor) + k / T = 0
    // Weak form for primary: (u_primary - k*u_neighbor, test)
    case Moose::Element:
      r = (_convectionCoeff * (_u[_qp] - _neighbor_value[_qp]) + ( _conductivity * _grad_u[_qp] * _normals[_qp])) * _test[_i][_qp];
      //r = (_u[_qp] - _kCoeff * _neighbor_value[_qp]) * _test[_i][_qp];
      break;

    // Neighbor residual: -(u_primary - k * u_neighbor, test_neighbor),
    // negative sign because the integration direction is opposite.
    case Moose::Neighbor:
      //r = - (_u[_qp] - _kCoeff * _neighbor_value[_qp]) * _test_neighbor[_i][_qp];
      r = (_convectionCoeff * (_u[_qp] - _neighbor_value[_qp]) + ( _neighborConductivity * _grad_neighbor_value[_qp] * _normals[_qp])) * - _test_neighbor[_i][_qp]; 
      break;
  }
  return r;
}

Real
InterfaceConvection::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * (_convectionCoeff * _phi[_j][_qp] + ( _grad_phi[_j][_qp] * _conductivity * _normals[_qp]));
      break;
    case Moose::NeighborNeighbor:
      jac = -_test_neighbor[_i][_qp] * (_convectionCoeff * - _phi_neighbor[_j][_qp] + (_grad_phi_neighbor[_j][_qp] * _neighborConductivity * _normals[_qp]));
      break;
    case Moose::NeighborElement:
      jac = -_test_neighbor[_i][_qp] * (_convectionCoeff * (_phi[_j][_qp]));//_phi[_j][_qp];
      break;
    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * (_convectionCoeff * (- _phi_neighbor[_j][_qp]));//-_kCoeff * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}
