#include "FHRApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
FHRApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Remove the runtime screen output message: LEGACY MODES ENABLED...
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

FHRApp::FHRApp(InputParameters parameters) : MooseApp(parameters)
{
  FHRApp::registerAll(_factory, _action_factory, _syntax);
}

FHRApp::~FHRApp() {}

void
FHRApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"FHRApp"});
  Registry::registerActionsTo(af, {"FHRApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
FHRApp::registerApps()
{
 registerApp(FHRApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
FHRApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  FHRApp::registerAll(f, af, s);
}
extern "C" void
FHRApp__registerApps()
{
  FHRApp::registerApps();
}
