//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
//*
//* Engy-5310: Computational Continuum Transport Phenomena
//* UMass Lowell, Nuclear Chemical Engineering
//* https://github.com/dpploy/engy-5310

#include "DiffusionFlux.h"

registerMooseObject("FuelPlateApp", DiffusionFlux);

InputParameters
DiffusionFlux::validParams()
{
  InputParameters params = VectorAuxKernel::validParams();
  params.addClassDescription("Computes the diffusion flux vector.");
  params.addRequiredCoupledVar("field", "The dependent variable field to be used.");
  params.addRequiredParam<Real>("diffCoeff", "Diffusion coefficient");

  return params;
}

DiffusionFlux::DiffusionFlux(const InputParameters & parameters): 
    VectorAuxKernel(parameters),
    _gradU(coupledGradient("field")),
    _diffCoeff(getParam<Real>("diffCoeff"))
{
}

RealVectorValue
DiffusionFlux::computeValue()
{
  return - _diffCoeff * _gradU[_qp];
}
