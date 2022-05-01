

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SourceTerm.h"

/**
 * All MOOSE based object classes you create must be registered using this macro. 
 * The first argument is the name of the App you entered in when running the stork.sh 
 * script with an "App" suffix. If you ran "stork.sh Example", then the argument here 
 * becomes "ExampleApp". The second argument is the name of the C++ class you created.
 
 Replace X
 */
registerMooseObject("Engy5310p2App", SourceTerm);

/**
 * This function defines the valid parameters for
 * this Kernel and their default values
 */
template<>
InputParameters validParams<SourceTerm>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("The equation term ($$), with the weak form of $$.");
  params.addParam<Real>("massArea",0.0,"Source term");
  params.addParam<Real>("henryC",0.0,"Source coupled term");
  params.addParam<Real>("massC",0.0,"Source term");
  params.addParam<Real>("heC",0.0,"Source coupled term");
  params.addRequiredCoupledVar("coupledVariable", "The variable field.");
  
  return params;
}

SourceTerm::SourceTerm(const InputParameters & parameters) : Kernel(parameters),
    // Set the coefficient for the equation term
    _massArea(getParam<Real>("massArea")),
	_henryC(getParam<Real>("henryC")),
	_massC(getParam<Real>("massC")),
	_heC(getParam<Real>("heC")),
	_uCoupled(coupledValue("coupledVariable"))
	
{
}

Real
SourceTerm::computeQpResidual()
{
 
return ((_massArea*(_u[_qp] - _henryC*_uCoupled[_qp])) - (_massC*(_uCoupled[_qp]-_heC*_u[_qp])))* _test[_i][_qp];

}

Real
SourceTerm::computeQpJacobian()
{
  
  return 0.0;
  
}
