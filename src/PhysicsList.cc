/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        PhysicsList.cc                                 /////
//      Description: An implementation of the Physics List class    /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////
#include "PhysicsList.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh" 
#include "G4EmStandardPhysics.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4EmProcessOptions.hh"
#include "PhysListMessenger.hh"

PhysicsList::PhysicsList(G4int ver):
  G4VModularPhysicsList(),
  messenger(new PhysListMessenger(this)),
  verbose(ver)
{
  defaultCutValue = 0.001*mm;//0.1*mm;
  AddFCPPhysics(0.01,100*MeV);//Default charge and mass of FCP
  RegisterPhysics( new G4EmStandardPhysics());
  RegisterPhysics( lPhys );
  SetVerboseLevel(verbose);
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetVerbose(verbose);
  param->SetMinEnergy(2e-3*eV);
  param->SetMaxEnergy(100*PeV);
}
PhysicsList::~PhysicsList()
{
}
void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
  //Construct FCP
  lPhys->ConstructParticle();
}


void PhysicsList::SetCuts()
{
  if (verboseLevel >0)
    {
    	G4cout << "PhysicsList::SetCuts:";
    	G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
    }
    
  //set cut values for gamma at first and for e- second and next for e+,
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(50*eV,100*PeV);
}
  
