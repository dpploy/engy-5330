//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InitialGuess.h"

registerMooseObject("SteamerApp", InitialGuess);

template <>
InputParameters
validParams<InitialGuess>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<Real>("coefficient", "The value of the initial condition");
  params.addParam<Real>("bias", "The value of the initial condition");
  return params;
}

InitialGuess::InitialGuess(const InputParameters & parameters): 
    InitialCondition(parameters), 
    _coefficient(getParam<Real>("coefficient")),
    _bias(getParam<Real>("bias"))
{
}

Real
InitialGuess::value(const Point & p)
{
  // The Point class is defined in libMesh. The spatial coordinates x,y,z can be accessed
  // individually using the parenthesis operator and a numeric index from 0..2
  // Linear wrt x initial guess
  return _bias + _coefficient * std::abs(p(0));
}
