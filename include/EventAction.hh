/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        RunAction.hh                                   /////
//      Description: Breif definition of the EventAction class      /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TRandom3.h"
#include <fstream>
using namespace std;

class RunAction;

/// Event action class
///

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction* runAction);
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);
  void AddEdep(G4double edep) { fEdep += edep; }
  void AddKEDiff(G4double kEDiff) { fdKE += kEDiff; }
  void Adddx(G4double dX) { fdx += dX; }

private:
  RunAction* fRunAction;
  G4int        EventID;
  G4double     Energy;
  G4double     Mass;
  G4double     Charge;
  G4double     fEdep;
  G4double     fdKE;
  G4double     fdx;
  TFile *f;
  TTree *t;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
