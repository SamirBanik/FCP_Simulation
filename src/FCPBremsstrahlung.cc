// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPBremsstrahlung.cc                                 //
//  Description: Radiation process for FCP particle.                  //
//                                                                    //
//  Adapted from G4MuBremsstrahlung by Samir Banik.                   //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
//  20180319  M. Kelsey -- Drop use of EmModel(i): broken in G4 10.4. //
////////////////////////////////////////////////////////////////////////

#include "FCPBremsstrahlung.hh"
#include "G4FCPMinus.hh"
#include "G4FCPPlus.hh"
#include "FCPBremsstrahlungModel.hh"
#include "G4EmParameters.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

FCPBremsstrahlung::FCPBremsstrahlung(const G4String& name)
  : G4VEnergyLossProcess(name),
    lowestKinEnergy(100.*eV),
    isInitialised(false)
{
  SetProcessSubType(fBremsstrahlung);
  SetSecondaryParticle(G4Gamma::Gamma());
  SetIonisation(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCPBremsstrahlung::~FCPBremsstrahlung()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCPBremsstrahlung::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4FCPPlus::Definition() || &p == G4FCPMinus::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double FCPBremsstrahlung::MinPrimaryEnergy(const G4ParticleDefinition*,
					      const G4Material*,
					      G4double) {
  return lowestKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPBremsstrahlung::
InitialiseEnergyLossProcess(const G4ParticleDefinition*,
			    const G4ParticleDefinition*) {
  if (isInitialised) return;		// Avoid unnecessary work

  isInitialised = true;
  G4VEmModel* theModel = new FCPBremsstrahlungModel;
  G4EmParameters* param = G4EmParameters::Instance();
  theModel->SetLowEnergyLimit(param->MinKinEnergy());
  theModel->SetHighEnergyLimit(param->MaxKinEnergy());
  AddEmModel(1, theModel, 0);		// No fluctuations
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPBremsstrahlung::PrintInfo()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

