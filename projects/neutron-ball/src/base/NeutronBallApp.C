#include "NeutronBallApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
NeutronBallApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Remove the runtime screen output message: LEGACY MODES ENABLED...
  params.set<bool>("use_legacy_material_output") = false;

  return params;
}

NeutronBallApp::NeutronBallApp(InputParameters parameters) : MooseApp(parameters)
{
  NeutronBallApp::registerAll(_factory, _action_factory, _syntax);
}

NeutronBallApp::~NeutronBallApp() {}

void
NeutronBallApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"NeutronBallApp"});
  Registry::registerActionsTo(af, {"NeutronBallApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
NeutronBallApp::registerApps()
{
  registerApp(NeutronBallApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
NeutronBallApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  NeutronBallApp::registerAll(f, af, s);
}
extern "C" void
NeutronBallApp__registerApps()
{
  NeutronBallApp::registerApps();
}
