//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArtificialDiffusion.h"

registerMooseObject("SteamerApp", ArtificialDiffusion);

template<>
InputParameters validParams<ArtificialDiffusion>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Artifical diffusion term.");
  params.addParam<Real>("diffCoeff",1.0,"Drift flux diffusion coefficient");
  return params;
}

ArtificialDiffusion::ArtificialDiffusion(const InputParameters & parameters): 
    Kernel(parameters),
    _diffCoeff(getParam<Real>("diffCoeff"))
{
}

Real ArtificialDiffusion::computeQpResidual()
{
 const RealVectorValue & vPrime = _grad_u[_qp];
 const RealVectorValue & thetaPrime = _grad_test[_i][_qp];

 return _diffCoeff * vPrime * thetaPrime;
}

Real
ArtificialDiffusion::computeQpJacobian()
{
 return 0.0;
}
