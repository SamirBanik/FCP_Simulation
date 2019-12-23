/////////////////////////////////////////////////////////////////////////
// ----------------------------------------------------------------------
//      File:        FCPPhyics.hh                                   /////
//      Description: Register process and models for FCPs           /////
//      21 December 2019, Samir Banik                               /////
// **********************************************************************
// ----------------------------------------------------------------------
/////////////////////////////////////////////////////////////////////////

#ifndef FCPPhysics_h
#define FCPPhysics_h 1

#include "G4VPhysicsConstructor.hh"
#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

class FCPPhysics : public G4VPhysicsConstructor
{
private:
  G4double lipcharge,lipmass;
public:
  FCPPhysics(G4double chg=0.01, G4double mass=100.*MeV);
  virtual ~FCPPhysics();
  virtual void ConstructParticle();
  virtual void ConstructProcess();
};
#endif
