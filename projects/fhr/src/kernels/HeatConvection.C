//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatConvection.h"
#include "Function.h"

registerMooseObject("FHRApp", HeatConvection);

template<>
InputParameters validParams<HeatConvection>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Heat convection term.");
  params.addParam<Real>("massDensity",1.0,"Fluid mass density");
  params.addRequiredParam<Real>("heatCapacity","Fluid heat capacity");
  params.addRequiredParam<FunctionName>("velocity","Fluid velocity Function");
  params.addRequiredParam<Real>("vmax", "Maximum Fluid Velocity");
  params.addRequiredParam<Real>("radius", "Channel Radius");
  params.addRequiredParam<Real>("n", "Inverse Velocity Profile Number");
  return params;
}

HeatConvection::HeatConvection(const InputParameters & parameters):
    Kernel(parameters),
    _massDensity(getParam<Real>("massDensity")),
    _heatCapacity(getParam<Real>("heatCapacity")),
    _func(getFunction("velocity")),
    _vmax(getParam<Real>("vmax")),
    _radius(getParam<Real>("radius")),
    _n(getParam<Real>("n"))
{
}

Real
HeatConvection::computeQpResidual()
{
  // calculate the velocity
  Real z_velocity = _vmax * pow((1 - (_func.value(_t, _q_point[_qp]) / _radius)), _n);
  // Put it in vector form
  RealVectorValue _vec = (0, 0, z_velocity);
  //std::cout<< _velocity_vec, '\n';
  //std::cout<< z_velocity;
  // Residual
  Real testing = _massDensity * _heatCapacity * _grad_u[_qp] * _vec * _test[_i][_qp];
  //std::cout<< testing;
  return testing;
}

Real
HeatConvection::computeQpJacobian()
{
  Real z_velocity = _vmax * pow((1 - (_func.value(_t, _q_point[_qp]) / _radius)), _n);
  RealVectorValue _vec = (0, 0, z_velocity);
  // Jacobian diagonal
  return _massDensity * _heatCapacity * _grad_phi[_j][_qp] * _vec * _test[_i][_qp];
}
