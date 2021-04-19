#include "SteamerApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
SteamerApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

SteamerApp::SteamerApp(InputParameters parameters) : MooseApp(parameters)
{
  SteamerApp::registerAll(_factory, _action_factory, _syntax);
}

SteamerApp::~SteamerApp() {}

void
SteamerApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"SteamerApp"});
  Registry::registerActionsTo(af, {"SteamerApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
SteamerApp::registerApps()
{
  registerApp(SteamerApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
SteamerApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SteamerApp::registerAll(f, af, s);
}
extern "C" void
SteamerApp__registerApps()
{
  SteamerApp::registerApps();
}
