//
// $Id: EventAction.cc December 23, 2019 $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4ParticleMomentum.hh"
#include <fstream>
#include "TMath.h"
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(),
    fRunAction(runAction),
    EventID(0),    
    Energy(0),
    Mass(0.),
    Charge(0.),
    fEdep(0.),
    fdKE(0.),
    fdx(0.)
{
  f = new TFile("DataFile.root","recreate");
  t = new TTree("tr","FCP Event data");
  //t->SetMaxTreeSize(1900000000);
  t->Branch("eID",&EventID,"eID/I");
  t->Branch("PrimaryMass",&Mass,"PrimaryMass/D");
  t->Branch("PrimaryCharge",&Charge,"PrimaryCharge/D");
  t->Branch("PrimaryKE",&Energy,"PrimaryKE/D");
  t->Branch("Eloss",&fdKE,"Eloss/D");
  t->Branch("Edep",&fEdep,"Edep/D");
  t->Branch("PathLength",&fdx,"PathLength/D");
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
  f->cd();
  t->Write();
  f->Close();
  //delete f;
}


void EventAction::BeginOfEventAction(const G4Event* evt)
{
  Energy = 0;
  fEdep = 0.;
  fdKE = 0.;
  fdx = 0.;
  EventID = evt->GetEventID();
  if(evt->GetEventID()%100000==0)
    {
      G4cout<<"Start of event # --> "<<evt->GetEventID()<<G4endl;
    }
  
  const PrimaryGeneratorAction* generatorAction
    = static_cast<const PrimaryGeneratorAction*>
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if (generatorAction)
    {
      const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
      Energy = particleGun->GetParticleEnergy()/MeV;
      Mass = particleGun->GetParticleDefinition()->GetPDGMass()/MeV;
      Charge = particleGun->GetParticleDefinition()->GetPDGCharge();
    }
}

//....oooOO0OOooo........oooOO0OOooo0........oooOO0OOooo........oooOO0OOooo......
void EventAction::EndOfEventAction(const G4Event*)
{   
  // accumulate statistics in run action
  fRunAction->AddEdep(fEdep);
  fdx       = fdx/cm;
  t->Fill();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
