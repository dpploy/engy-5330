#include "HeatedPlateApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
HeatedPlateApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

HeatedPlateApp::HeatedPlateApp(InputParameters parameters) : MooseApp(parameters)
{
  HeatedPlateApp::registerAll(_factory, _action_factory, _syntax);
}

HeatedPlateApp::~HeatedPlateApp() {}

void
HeatedPlateApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"HeatedPlateApp"});
  Registry::registerActionsTo(af, {"HeatedPlateApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
HeatedPlateApp::registerApps()
{
  registerApp(HeatedPlateApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
HeatedPlateApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  HeatedPlateApp::registerAll(f, af, s);
}
extern "C" void
HeatedPlateApp__registerApps()
{
  HeatedPlateApp::registerApps();
}
