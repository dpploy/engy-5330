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

#include "InterfaceFlux.h"

registerMooseObject("FHRApp", InterfaceFlux);

InputParameters
InterfaceFlux::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addParam<Real>("conductivity",1.0,"Primary diffusion coefficient.");
  params.addParam<Real>("neighborConductivity",1.0,"Neighboring diffusion coefficient.");
  return params;
}

InterfaceFlux::InterfaceFlux(const InputParameters & parameters):
    InterfaceKernel(parameters),
    _conductivity(getParam<Real>("conductivity")),
    _neighborConductivity(getParam<Real>("neighborConductivity"))
{
}

Real
InterfaceFlux::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  // While on the primary side compute quantities on the neighbor side
  RealVectorValue qNeighbor = - _neighborConductivity * _grad_neighbor_value[_qp];
  RealVectorValue normalNeighbor = - _normals[_qp];
  Real qNeighborN = qNeighbor * normalNeighbor;

  // While on the neighbor side compute quantities on the primary side
  RealVectorValue qElement = - _conductivity * _grad_u[_qp];
  RealVectorValue normal = _normals[_qp];
  Real qElementN = qElement * normal;

  switch (type)
  {
    case Moose::Element:

      r = - qNeighborN * _test[_i][_qp];
      break;

    case Moose::Neighbor:

      r = - qElementN * _test_neighbor[_i][_qp];
      break;
  }

  return r;
}

Real
InterfaceFlux::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
    case Moose::NeighborNeighbor:
      break;

    case Moose::NeighborElement:
      jac = _conductivity * _grad_phi[_j][_qp] * _normals[_qp] * _test_neighbor[_i][_qp]; 
      break;

    case Moose::ElementNeighbor:
      jac = - _neighborConductivity * _grad_phi_neighbor[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      break;
  }

  return jac;
}
