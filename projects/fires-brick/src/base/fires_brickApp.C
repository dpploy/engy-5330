#include "fires_brickApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
fires_brickApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

fires_brickApp::fires_brickApp(InputParameters parameters) : MooseApp(parameters)
{
  fires_brickApp::registerAll(_factory, _action_factory, _syntax);
}

fires_brickApp::~fires_brickApp() {}

void
fires_brickApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"fires_brickApp"});
  Registry::registerActionsTo(af, {"fires_brickApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
fires_brickApp::registerApps()
{
  registerApp(fires_brickApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
fires_brickApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  fires_brickApp::registerAll(f, af, s);
}
extern "C" void
fires_brickApp__registerApps()
{
  fires_brickApp::registerApps();
}
