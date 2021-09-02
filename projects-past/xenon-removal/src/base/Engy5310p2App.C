#include "Engy5310p2App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
Engy5310p2App::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

Engy5310p2App::Engy5310p2App(InputParameters parameters) : MooseApp(parameters)
{
  Engy5310p2App::registerAll(_factory, _action_factory, _syntax);
}

Engy5310p2App::~Engy5310p2App() {}

void
Engy5310p2App::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"Engy5310p2App"});
  Registry::registerActionsTo(af, {"Engy5310p2App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
Engy5310p2App::registerApps()
{
  registerApp(Engy5310p2App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
Engy5310p2App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Engy5310p2App::registerAll(f, af, s);
}
extern "C" void
Engy5310p2App__registerApps()
{
  Engy5310p2App::registerApps();
}
