// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPIonisation.cc                                     //
//  Description: Ionisation and delta-ray production for LIP particle //
//                                                                    //
//  Adapted from G4MuIonisation by Samir Banik.                       //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
//  20180319  M. Kelsey -- Drop use of EmModel(i): broken in G4 10.4. //
////////////////////////////////////////////////////////////////////////

#include "FCPIonisation.hh"
#include "G4FCPMinus.hh"
#include "G4FCPPlus.hh"
#include "FCPBetheBlochModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4BohrFluctuations.hh"
#include "G4BraggModel.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4ICRU73QOModel.hh"
#include "G4IonFluctuations.hh"
#include "G4PAIModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UniversalFluctuation.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

FCPIonisation::FCPIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(0),
    theBaseParticle(0),
    isInitialised(false)
{
  mass = ratio = 0;
  SetProcessSubType(fIonisation);
  SetSecondaryParticle(G4Electron::Electron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCPIonisation::~FCPIonisation() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCPIonisation::IsApplicable(const G4ParticleDefinition& p) {
  return (&p == G4FCPPlus::Definition() || &p == G4FCPMinus::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double FCPIonisation::MinPrimaryEnergy(const G4ParticleDefinition*,
					  const G4Material*,
					  G4double cut) {
  G4double x = 0.5*cut/electron_mass_c2;
  G4double gam = x*ratio + std::sqrt((1. + x)*(1. + x*ratio*ratio));
  return mass*(gam - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPIonisation::
InitialiseEnergyLossProcess(const G4ParticleDefinition* part,
			    const G4ParticleDefinition* bpart) {
  if (isInitialised) return;		// Avoid unnecessary work

  theParticle = part;
  theBaseParticle = bpart;
  
  mass = theParticle->GetPDGMass();
  ratio = electron_mass_c2/mass;
  
  G4EmParameters* param = G4EmParameters::Instance();
  G4double elow = param->MinKinEnergy();//0.2*MeV;
  G4double emax = param->MaxKinEnergy();//param->MaxKinEnergy();
  //G4double ehigh = std::min(4e-3*eV, emax);//std::min(1*GeV, emax);

  G4PAIModel* theModel = new G4PAIModel(theParticle,"PAIModel");
  //G4PAIModel* theFluct = new G4PAIModel(theParticle,"PAIModel");
  //%%%%%%%%%%%%%%%% Region for PAIModel added July 26%%
  const G4RegionStore* theRegionStore = G4RegionStore::GetInstance();
  G4Region* Region = theRegionStore->GetRegion("DefaultRegionForTheWorld");
  //theModel = new G4PAIModel(theParticle,"PAIModel");
  //theFluct = new G4PAIModel(theParticle,"PAIModel");//G4UniversalFluctuation;
  theModel->SetLowEnergyLimit(elow);
  theModel->SetHighEnergyLimit(emax);
  //theFluct->SetLowEnergyLimit(elow);
  //theFluct->SetHighEnergyLimit(emax);
  //FCP_PAIModel->SetVerboseLevel(0);
  AddEmModel(0,theModel,theModel,Region);
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPIonisation::PrintInfo() {;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




