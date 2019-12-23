/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        SteppingAction.cc                              /////
//      Description: An implementation of the SteppingAction class  /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4VITProcess.hh"
#include <fstream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringFiducial(0),
  fWorldPhysical(0),
  volume1(0),
  volume2(0),
  fParticleGun(0),
  PrimaryPDG(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringFiducial){
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
      (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringFiducial = detectorConstruction->GetScoringFiducial();
  }
  // if (!fWorldPhysical) {
  //   fWorldPhysical = detectorConstruction->GetWorldPhysical();
  // }
  
  if(!fParticleGun){
    const PrimaryGeneratorAction* generatorAction
      = static_cast<const PrimaryGeneratorAction*>
      (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    fParticleGun = generatorAction->GetParticleGun();
    PrimaryPDG = fParticleGun->GetParticleDefinition()->GetPDGEncoding();
  }
  G4StepPoint* point1 = step->GetPreStepPoint();
  G4StepPoint* point2 = step->GetPostStepPoint();
  // get volume of the current step
  volume1 = point1->GetPhysicalVolume();
  //volume2 = point2->GetPhysicalVolume();
  if(volume1  != fScoringFiducial) return;

  
  G4Track* track         = step->GetTrack();
  G4int PDG    = track->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();
  G4double edepStep = step->GetTotalEnergyDeposit();
  G4double StepSize = track->GetStepLength();
  G4double PreStepKE = point1->GetKineticEnergy();
  G4double PostStepKE = point2->GetKineticEnergy();
  G4double KEDiff = PreStepKE - PostStepKE;

  // collect energy deposited in this step
  fEventAction->AddEdep(edepStep);

  //  collect energy lost in this step by Primary
  if(PDG == PrimaryPDG)
    {  
      fEventAction->AddKEDiff(KEDiff);
      fEventAction->Adddx(StepSize);
    }
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

