// $Id: 
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        FCPMultipleScattering.cc                             //
//                                                                    //
//  Adapted from G4MuMultipleScattering by Samir Banik.               //
//                                                                    //
//  Author:      Samir Banik                                          //
//  Date:        14-th June 2018                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "FCPMultipleScattering.hh"
#include "G4SystemOfUnits.hh"
#include "G4WentzelVIModel.hh"
#include "G4MscStepLimitType.hh"
#include "G4FCPPlus.hh"
#include "G4FCPMinus.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

FCPMultipleScattering::FCPMultipleScattering(const G4String& pnam)
  : G4VMultipleScattering(pnam)
{
  isInitialized = false;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FCPMultipleScattering::~FCPMultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool FCPMultipleScattering::IsApplicable (const G4ParticleDefinition& p)
{
  //return (p.GetPDGCharge() != 0.0 && !p.IsShortLived());
  return (&p == G4FCPPlus::Definition() || &p == G4FCPMinus::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCPMultipleScattering::InitialiseProcess(const G4ParticleDefinition*)
{
  // Modification of parameters between runs
  if(isInitialized) { return; }
  /*
  if(!EmModel(1)) { SetEmModel(new G4WentzelVIModel(), 1); } //Samir 07.05.18
  AddEmModel(1, EmModel(1));
  isInitialized = true;
  */
  G4VEmModel* theModel = 0;
  theModel = new G4WentzelVIModel();
  AddEmModel(0, theModel);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCPMultipleScattering::PrintInfo()
{
  G4cout << "      RangeFactor= " << RangeFactor()
         << ", step limit type: " << StepLimitType()
         << ", lateralDisplacement: " << LateralDisplasmentFlag()
	 << ", polarAngleLimit(deg)= " << PolarAngleLimit()/degree
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

