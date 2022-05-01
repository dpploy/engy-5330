//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArtificialViscosity.h"

registerMooseObject("SteamerApp", ArtificialViscosity);

template<>
InputParameters validParams<ArtificialViscosity>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Artifical viscosity term.");
  params.addParam<Real>("viscosity",1.0,"Artificial viscosity coefficient");
  return params;
}

ArtificialViscosity::ArtificialViscosity(const InputParameters & parameters): 
    Kernel(parameters),
    _viscosity(getParam<Real>("viscosity"))
{
}

Real ArtificialViscosity::computeQpResidual()
{
 const RealVectorValue & vPrime = _grad_u[_qp];
 const RealVectorValue & thetaPrime = _grad_test[_i][_qp];

 return _viscosity * vPrime * thetaPrime;
}

Real
ArtificialViscosity::computeQpJacobian()
{
 return 0.0;
}
