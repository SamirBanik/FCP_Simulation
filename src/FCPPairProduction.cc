// $Id$
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPPairProduction.cc                                 //
//  Description: e+e- pair production by FCPs                         //
//                                                                    //
//  Adapted from G4MuPairProduction by Samir Banik.                   //
//                                                                    //
//  Author:	 Samir Banik                                          //
//  Date:        9-th May 2017                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////
#include "FCPPairProduction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Positron.hh"
#include "G4VEmModel.hh"
#include "FCPPairProductionModel.hh"
#include "G4ElementData.hh"
#include "G4EmParameters.hh"
#include "G4FCPPlus.hh"
#include "G4FCPMinus.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

FCPPairProduction::FCPPairProduction(const G4String& name)
  : G4VEnergyLossProcess(name),
    theParticle(0),
    lowestKinEnergy(1.*GeV),
    isInitialised(false)
{
  SetProcessSubType(fPairProdByCharged);
  SetSecondaryParticle(G4Positron::Positron());
  SetIonisation(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCPPairProduction::~FCPPairProduction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool FCPPairProduction::IsApplicable(const G4ParticleDefinition& p)
{
  //return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 10.0*MeV);
  return (&p == G4FCPPlus::Definition() || &p == G4FCPMinus::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double FCPPairProduction::MinPrimaryEnergy(const G4ParticleDefinition*,
					      const G4Material*,
					      G4double)
{
  return lowestKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPPairProduction::InitialiseEnergyLossProcess(
                         const G4ParticleDefinition* part,
			 const G4ParticleDefinition*)
{
  if (!isInitialised) {
    isInitialised = true;

    theParticle = part;

    if (!EmModel()) { SetEmModel(new FCPPairProductionModel(part)); }

    G4double limit = part->GetPDGMass()*8;
    if(limit > lowestKinEnergy) { lowestKinEnergy = limit; }

    G4VEmFluctuationModel* fm = 0;
    G4EmParameters* param = G4EmParameters::Instance();
    EmModel()->SetLowEnergyLimit(param->MinKinEnergy());
    EmModel()->SetHighEnergyLimit(param->MaxKinEnergy());
    AddEmModel(1, EmModel(), fm);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCPPairProduction::PrintInfo()
{
  G4ElementData* ed = EmModel()->GetElementData();
  if(ed) {
    for(G4int Z=1; Z<93; ++Z) {
      G4Physics2DVector* pv = ed->GetElement2DData(Z);
      if(pv) {
        G4cout << "      Sampling table " << pv->GetLengthY()
	       << "x" << pv->GetLengthX() << "; from "
	       << exp(pv->GetY(0))/GeV << " GeV to " 
	       << exp(pv->GetY(pv->GetLengthY()-1))/TeV 
	       << " TeV " << G4endl;
	break;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




