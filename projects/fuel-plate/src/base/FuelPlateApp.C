#include "FuelPlateApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
FuelPlateApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

FuelPlateApp::FuelPlateApp(InputParameters parameters) : MooseApp(parameters)
{
  FuelPlateApp::registerAll(_factory, _action_factory, _syntax);
}

FuelPlateApp::~FuelPlateApp() {}

void
FuelPlateApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"FuelPlateApp"});
  Registry::registerActionsTo(af, {"FuelPlateApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
FuelPlateApp::registerApps()
{
  registerApp(FuelPlateApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
FuelPlateApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FuelPlateApp::registerAll(f, af, s);
}
extern "C" void
FuelPlateApp__registerApps()
{
  FuelPlateApp::registerApps();
}
