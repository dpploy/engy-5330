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

#include "InterfaceJump.h"

registerMooseObject("HeatedPlateApp", InterfaceJump);

defineLegacyParams(InterfaceJump);

InputParameters
InterfaceJump::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addClassDescription("Wall heat convection transport");
  params.addRequiredParam<Real>("transferCoeff", "Heat transfer coefficient");
  params.addRequiredParam<Real>("adsorptionCoeff", "Adsorption coefficient");
  params.addRequiredParam<Real>("thermCondCoeff", "Primary thermal conductivity coefficient");
  return params;
}

InterfaceJump::InterfaceJump(const InputParameters & parameters): 
  InterfaceKernel(parameters), 
  _transferCoeff(getParam<Real>("transferCoeff")),
  _adsorptionCoeff(getParam<Real>("adsorptionCoeff")),
  _thermCondCoeff(getParam<Real>("thermCondCoeff"))
{
}

Real
InterfaceJump::computeQpResidual(Moose::DGResidualType type)
{
  Real C = _adsorptionCoeff;
  Real h = _transferCoeff;
  Real k = _thermCondCoeff;

  RealVectorValue q_primary = - k * _grad_u[_qp];
  RealVectorValue normal = _normals[_qp];
  Real qn_primary = q_primary * normal;

  Real r = 0;
  switch (type)
  {
    // Primary residual = (1-C) * u_neighbor - u_primary + 1/h * qn_primary
    // Weak form for primary: ( (1-C) * u_neighbor - u_primary + 1/h * qn_primary, test)
    case Moose::Element:

      r = - ((1.0-C) * _neighbor_value[_qp] - _u[_qp] + 1.0/h * qn_primary) * _test[_i][_qp];
      break;

    // Neighbor residual: - ((1-C) * u_neighbor - u_primary + 1/h * qn_primary, test_neighbor)
    // negative sign because the integration direction is opposite.
    case Moose::Neighbor:

      r = ((1.0-C) * _neighbor_value[_qp] - _u[_qp] + 1.0/h * qn_primary) * _test_neighbor[_i][_qp];
      break;
  }
  return r;
}

Real
InterfaceJump::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;
  switch (type)
  {
    case Moose::ElementElement:
      jac = _test[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::NeighborNeighbor:
      jac = -_test_neighbor[_i][_qp] * -1 * _phi_neighbor[_j][_qp];
      break;
    case Moose::NeighborElement:
      jac = -_test_neighbor[_i][_qp] * _phi[_j][_qp];
      break;
    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * 1 * _phi_neighbor[_j][_qp];
      break;
  }
  return jac;
}
