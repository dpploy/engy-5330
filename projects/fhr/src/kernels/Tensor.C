//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Tensor.h"

registerMooseObject("FHRApp", Tensor);

template<>
InputParameters validParams<Tensor>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($...$), with the weak form of $...$.");
  params.addParam<Real>("A",1.0,"Equation Term Coefficient");
  params.addParam<Real>("B");
  params.addParam<Real>("C");
  params.addParam<Real>("density");
  return params;
}

Tensor::Tensor(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    A(getParam<Real>("A")),
    B(getParam<Real>("B")),
    C(getParam<Real>("C")),
    density(getParam<Real>("density"));
{
}
Real gradTensor()
{
	return (4 * A^2 * x^3 + 6 * A * B * x^2 + 4 * A * C * x + 2 * B^2 * x + 2 * B * C);
}
Real
Tensor::computeQpResidual()
{
  // Implement the return
  tens = gradTensor();
  Real blin =  - density * gradTensor  *  _test[_i][_qp];
  //printf("%f ", _qp);
  return blin;
  //FIXME return 0.0 // remove this line
}

Real
Tensor::computeQpJacobian()
{
  // Implement the return
  //Real din = - diffCoeff * _grad_phi[_j][_qp] * _grad_test[_i][_qp];
  //printf("%f ", _qp);
  //return din;

  return 0.0
}
