/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        SteppingAction.hh                              /////
//      Description: Breif definition of the SteppingAction class   /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <fstream>
using namespace std;

class EventAction;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4ParticleGun;

/// Stepping action class
class SteppingAction : public G4UserSteppingAction
{
 public:
  SteppingAction(EventAction* eventAction);
  virtual ~SteppingAction();

  // method from the base class
  virtual void UserSteppingAction(const G4Step*);
  
 private:
  EventAction*  fEventAction;
  G4VPhysicalVolume* fScoringFiducial;
  G4VPhysicalVolume* fWorldPhysical;
  G4VPhysicalVolume* volume1;
  G4VPhysicalVolume* volume2;
  const G4ParticleGun* fParticleGun;

  G4int PrimaryPDG;
  ofstream fstep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
