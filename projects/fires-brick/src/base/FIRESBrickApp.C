#include "FIRESBrickApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
FIRESBrickApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy material output, i.e., output properties on INITIAL as well as TIMESTEP_END
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

FIRESBrickApp::FIRESBrickApp(InputParameters parameters) : MooseApp(parameters)
{
  FIRESBrickApp::registerAll(_factory, _action_factory, _syntax);
}

FIRESBrickApp::~FIRESBrickApp() {}

void
FIRESBrickApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAll(f, af, syntax);
  Registry::registerObjectsTo(f, {"FIRESBrickApp"});
  Registry::registerActionsTo(af, {"FIRESBrickApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
FIRESBrickApp::registerApps()
{
  registerApp(FIRESBrickApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
FIRESBrickApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FIRESBrickApp::registerAll(f, af, s);
}
extern "C" void
FIRESBrickApp__registerApps()
{
  FIRESBrickApp::registerApps();
}
