////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        PhysListMessenger.cc                                 //
//  Description: A PhysicsList messenger to specify Mass and Charge   //
//               of FCPs                                              //
//  Author:      Varchaswi K. S. Kashyap                              //
//  Date:        December 21, 2019                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "PhysListMessenger.hh"
#include "G4UIcommand.hh"
#include "PhysicsList.hh"
#include <sstream>
#include "G4UIcmdWithAString.hh"
PhysListMessenger::PhysListMessenger(PhysicsList* PhyList)
  : thePhysicsList(PhyList), LIPCmd(0)
{
  LIPCmd = new G4UIcommand("/FCP/Physics/ParticleProperties",this);
  LIPCmd->SetGuidance("Setup particle properties");

  //User command for charge
  G4UIparameter* q = new G4UIparameter("charge",'d',false);
  q->SetGuidance("Electric charge");
  q->SetDefaultValue("0.01");
  LIPCmd->SetParameter(q);


  //User command for mass
  G4UIparameter* mass = new G4UIparameter("mass",'d',false);
  mass->SetGuidance("mass");
  mass->SetParameterRange("mass>0.");
  mass->SetDefaultValue("100");
  LIPCmd->SetParameter(mass);

  G4UIparameter* unit = new G4UIparameter("unit",'s',false);
  LIPCmd->SetParameter(unit);
  unit->SetDefaultValue("MeV");
  LIPCmd->AvailableForStates(G4State_PreInit);
}

PhysListMessenger::~PhysListMessenger() {
  delete LIPCmd;
}

void PhysListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if (command == LIPCmd)
   { G4double q, mass;
     G4String unts;
     std::istringstream is(newValue);
     is >> q >> mass >> unts;
     G4String unit = unts;
     G4double vUnit = G4UIcommand::ValueOf(unit);
     thePhysicsList->AddFCPPhysics(q,mass*vUnit);
   }
}
