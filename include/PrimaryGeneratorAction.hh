///////////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------------
//      File:        PrimaryGeneratorAction.hh                            /////
//      Description: Breif definition of the PrimaryGeneratorAction class /////
//      21 December 2019, Samir Banik                                     /////
// ****************************************************************************
// ---------------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include <fstream>
#include "Randomize.hh"
using namespace std;

class G4ParticleGun;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 225.73 MeV gamma, randomly distribued 

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  virtual ~PrimaryGeneratorAction();
  
  // method from the base class
  virtual void GeneratePrimaries(G4Event*);         
  
  // method to access particle gun
  const G4ParticleGun* GetParticleGun() const { return fParticleGun; }
  

  
private:
  G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
